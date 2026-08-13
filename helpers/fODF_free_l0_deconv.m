function H = fODF_free_l0_deconv()
% H = fODF_free_l0_deconv()
%
% EXPERIMENTAL. One regularized spherical deconvolution that can be run in
% either of two conventions, so the only difference between them is the one
% being studied:
%
%   'fixed'  p_00 = 1 is imposed, and only l >= 2 is estimated. This is what
%            SMI.get_plm_from_S_and_kernel does, and it is what makes an SMI
%            fODF integrate to 1 in every voxel.
%   'free'   every coefficient including l = 0 is estimated. The fODF no
%            longer has unit mass and its amplitude carries density, as an
%            MRtrix FOD does.
%
%   H.solve(A, y, Ylm, reg, mode)   one voxel
%   H.reg_defaults(Lmax, reg)       fill in defaults and precompute Gamma
%
% A     [Ndwi x Nlm]  full design matrix, K_l(b) * Y_lm(u) * sqrt((2l+1)*4pi),
%                     INCLUDING the l = 0 column
% y     [Ndwi x 1]    the normalized signal, with NOTHING subtracted
% Ylm   [Ndirs x Nlm] fODF amplitude per unit plm, INCLUDING the l = 0 column
% reg   options, same fields as SMI.fODF_RegularizationDefaults
% mode  'fixed' or 'free'
%
% Returns p as the FULL [Nlm x 1] coefficient vector in both modes, with
% p(1) = 1 in 'fixed' mode, so the caller does not branch on the convention.
%
% WHY THIS EXISTS RATHER THAN A FLAG IN SMI.m. The two conventions differ in
% four places, not one, and putting them side by side in one function is the
% only way to be sure nothing else drifted between them:
%
%  1. the design matrix keeps or drops its l = 0 column;
%  2. the data has the l = 0 prediction subtracted, or does not;
%  3. the non-negativity constraint is INHOMOGENEOUS in 'fixed' mode -- the
%     fixed isotropic floor 1/(4pi) moves to the right-hand side -- and
%     HOMOGENEOUS in 'free' mode, where there is no floor at all. That is the
%     substantive difference: with p_00 fixed the fODF cannot be dimmed below
%     its floor, and with p_00 free it can go to zero;
%  4. the amplitude threshold tau*mean(fODF) is a constant in 'fixed' mode and
%     depends on the current solution in 'free' mode.
%
% VALIDATION. In 'fixed' mode this function must reproduce
% SMI.get_plm_from_S_and_kernel to machine precision. Anything that uses the
% 'free' results without checking that first is trusting a reimplementation,
% which this repository has been burned by before.
%
% Lives in its own file rather than as a local function because MATLAB requires
% a script's local functions at the END of the file while Octave cannot call
% them there at all. Same reasoning as helpers/fODF_peak_score.m.
H = struct();
H.solve        = @solve_voxel;
H.reg_defaults = @reg_defaults;
end

% =====================================================================
function reg = reg_defaults(Lmax, reg)
% Mirrors SMI.fODF_RegularizationDefaults, with one deliberate difference:
% Gamma covers ALL coefficients including l = 0, and its l = 0 entry is always
% zero. Tikhonov damping exists to suppress the high orders the kernel
% attenuates most; shrinking p_00 would shrink the density estimate itself,
% which is the quantity 'free' mode exists to measure.
if nargin < 2 || isempty(reg), reg = struct(); end
d = struct('flag_nonneg',0,'lambda_nonneg',1,'tau',0.1,'Niter',50, ...
           'Ndirs',300,'Lmax_init',4,'max_neg_fraction',0.9, ...
           'lambda_tikhonov',0,'TikhonovMatrix','identity');
fn = fieldnames(d);
for ii = 1:numel(fn)
    if ~isfield(reg,fn{ii}) || isempty(reg.(fn{ii})), reg.(fn{ii}) = d.(fn{ii}); end
end
L_all = repelem(0:2:Lmax, 2*(0:2:Lmax)+1);
switch lower(reg.TikhonovMatrix)
    case 'identity'
        g = ones(size(L_all));
    case {'laplacebeltrami','laplace-beltrami','lb'}
        g = L_all.*(L_all+1);
        g = g/max(g);
    otherwise
        error('TikhonovMatrix must be ''identity'' or ''laplacebeltrami''');
end
g(L_all == 0) = 0;                       % never damp the density term
reg.Gamma_full = diag(g);
reg.L_all      = L_all;
reg.init_full  = L_all <= reg.Lmax_init; % includes l = 0
end

% =====================================================================
function [p, info] = solve_voxel(A, y, Ylm, reg, mode)
% info = [Niterations; Ndirections constrained; converged flag]
Nlm  = size(A,2);
free = strcmpi(mode,'free');
info = [0;0;1];
p    = zeros(Nlm,1);

if ~all(isfinite(A(:))) || ~all(isfinite(y))
    p = nan(Nlm,1); info(3) = 0; return
end

if free
    Ad  = A;                 % keep the l = 0 column
    yd  = y(:);
    Yd  = Ylm;               % amplitude includes the l = 0 contribution
    c0  = zeros(size(Ylm,1),1);   % no fixed floor: the constraint is homogeneous
    idx = 1:Nlm;
else
    Ad  = A(:,2:end);        % SMI's arrangement
    yd  = y(:) - A(:,1);     % p_00 = 1 moved to the right-hand side
    Yd  = Ylm(:,2:end);
    c0  = Ylm(:,1);          % the fixed isotropic floor, 1/(4*pi)
    idx = 2:Nlm;
end
Nc = size(Ad,2);

scale_A = norm(Ad,'fro')/sqrt(size(Ad,1));
if ~isfinite(scale_A) || scale_A < eps
    if ~free, p(1) = 1; end
    info(3) = 0; return
end

G = reg.Gamma_full(idx, idx);
if reg.lambda_tikhonov > 0
    Tik = (reg.lambda_tikhonov*scale_A)*G;
else
    Tik = zeros(0, Nc);
end
rhs_tik = zeros(size(Tik,1),1);

% ---- initial solution
if reg.flag_nonneg
    im = reg.init_full(idx);
else
    im = true(1,Nc);
end
q = zeros(Nc,1);
if size(Tik,1) > 0, Ti = Tik(im,im); else, Ti = zeros(0,sum(im)); end
q(im) = [Ad(:,im); Ti] \ [yd; zeros(size(Ti,1),1)];

if reg.flag_nonneg
    scale_Y = norm(Yd,'fro')/sqrt(size(Yd,1));
    w       = reg.lambda_nonneg*scale_A/scale_Y;
    Ndirs   = size(Yd,1);
    % In 'free' mode the mean fODF is not a constant, so the threshold is
    % recomputed from the current solution rather than fixed once.
    thr      = reg.tau*mean(c0 + Yd*q);
    neg_prev = [];
    conv     = 0;
    neg      = false(Ndirs,1);
    for it = 1:reg.Niter
        neg = (c0 + Yd*q) < thr;
        if ~any(neg) || (it > 1 && isequal(neg, neg_prev)), conv = 1; break, end
        if sum(neg) > reg.max_neg_fraction*Ndirs, break, end
        L = w*Yd(neg,:);
        q = [Ad; L; Tik] \ [yd; -w*c0(neg); rhs_tik];
        if free
            thr = reg.tau*mean(Yd*q);      % track the density as it changes
        end
        neg_prev = neg; info(1) = it;
    end
    info(2) = sum(neg); info(3) = conv;
end

if free
    p = q(:);
else
    p = [1; q(:)];
end
end
