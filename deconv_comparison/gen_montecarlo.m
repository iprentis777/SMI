function gen_montecarlo(SNR, NREP, tag)
% gen_montecarlo(SNR, NREP, tag)
%
% The Monte Carlo experiment of REPORT_SMI_deconvolution_MonteCarlo.md, in the
% design of Jeurissen et al. (2014): synthesise noise-free signal vectors for
% crossing white matter fibres by forward convolution, add complex Gaussian
% noise so the magnitude is Rician, and deconvolve NREP independent
% realisations of each condition.
%
% Conditions: a single fibre plus 15, 45 and 60 degree crossings of two equal
% fibre populations. Each fibre population is a Watson distribution, kappa =
% KAPPA, so the ground truth carries realistic dispersion rather than being a
% delta -- which matters, because a response function estimated from real white
% matter absorbs dispersion and a delta ground truth would hand every method a
% mismatch that does not exist in practice.
%
% Written per SNR:
%   mc_dwi_<tag>        [GRID x Ndwi]  the noisy signal, byte identical for
%                                      every method
%   sh_smiC_<tag>       [NVOX x Nlm]   SMI fODF, constrained deconvolution
%   sh_smiU_<tag>       [NVOX x Nlm]   SMI fODF, unconstrained (plain LLS)
%   kernel_<tag>        [NVOX x Nk]    the SMI kernel, shared by both
%   reginfo_<tag>       [NVOX x 3]     iterations / constrained dirs / converged
%
% The two SMI arms share ONE kernel fit: SMI.fit is run once and the
% unconstrained deconvolution is then recomputed from out.kernel by calling
% SMI.get_plm_from_S_and_kernel directly. Any difference between the arms is
% therefore the deconvolution alone, not a different kernel.
%
% Coefficients are stored as SH coefficients in SMI's convention
% (f_lm = [1 plm].*sqrt((2l+1)/(4pi)), so the fODF integrates to 1) together
% with the SH matrix Y_smi that generated them, so the Python side can map them
% into dipy's basis exactly rather than by assuming the two agree.

more off
IO = binio();

bvals     = IO.load('bvals'); bvals = bvals(:)';
bvecs     = IO.load('bvecs');
eval_dirs = IO.load('eval_dirs');
Ndwi      = numel(bvals);

LMAX_FIT = 6;      % every method deconvolves at this angular order
LMAX_GT  = 8;      % the truth is NOT band limited to it (8 is the SM kernel's
                   % ceiling, SMI.m:2470-2480)
CS       = 1;
D_FW     = 3;
KAPPA    = 16;
ANGLES   = [0 15 45 60];        % 0 = single fibre
K_WM     = [0.60 2.0 2.0 0.50 0.02];

H  = fODF_modulation_helpers();
dq = H.dirs(3000);

n1 = [0.30 -0.50 0.81]; n1 = n1/norm(n1);

NCOND = numel(ANGLES);
NVOX  = NCOND*NREP;
GRID  = pick_grid(NVOX);
fprintf('SNR %g, %d conditions x %d reps = %d voxels, grid %s\n', ...
        SNR, NCOND, NREP, NVOX, mat2str(GRID));

% --------------------------------------------------- ground truth signals
Y_gt  = SMI.get_even_SH(eval_dirs, LMAX_GT, CS);
L_gt  = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';
keep6 = L_gt <= LMAX_FIT;

S_cond  = zeros(NCOND, Ndwi);
sh_gt   = zeros(NCOND, numel(L_gt));           % truth, full order
sh_gt6  = zeros(NCOND, sum(keep6));            % truth band limited to the fit
ax      = nan(2, 3, NCOND);
for ic = 1:NCOND
    if ANGLES(ic) == 0
        axes_ = {n1};
    else
        axes_ = {n1, rotate_about(n1, ANGLES(ic))};
    end
    fod = zeros(size(dq,1),1);
    for k = 1:numel(axes_), fod = fod + H.watson(dq, axes_{k}, KAPPA); end
    plm_gt = H.mixture_plm(fod, dq, LMAX_GT, CS);
    coef   = [1; plm_gt(:)].*sqrt((2*L_gt+1)/(4*pi));
    sh_gt(ic,:)  = coef';
    sh_gt6(ic,:) = coef(keep6)';
    s = H.signal(plm_gt, [K_WM 1 1], bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                 bvecs, LMAX_GT, CS, D_FW);
    S_cond(ic,:) = s(:)';
    for k = 1:numel(axes_), ax(k,:,ic) = axes_{k}; end
end

cond_id = repelem((1:NCOND)', NREP, 1);
S_clean = S_cond(cond_id, :);

% ------------------------------------------------------------ Rician noise
% S0 = 1 by construction, so sigma = 1/SNR. Complex Gaussian noise on a real
% signal, magnitude taken: exactly Rician.
sigma = 1/SNR;
rand('seed', 31415); randn('seed', 31415);
S_noisy = sqrt((S_clean + sigma*randn(size(S_clean))).^2 + ...
               (         sigma*randn(size(S_clean))).^2);
dwi = reshape(S_noisy, [GRID Ndwi]);

% Everything is tagged, including the parts that do not depend on the noise
% level: the SNR arms are run as concurrent processes and two of them writing
% the same file is a race, not a saving.
IO.save(['mc_dwi_' tag], dwi);
IO.save(['mc_cond_id_' tag], cond_id);
IO.save(['mc_gt_axes_' tag], ax);
IO.save(['mc_sh_gt_' tag], sh_gt);
IO.save(['mc_sh_gt6_' tag], sh_gt6);
IO.save(['mc_angles_' tag], ANGLES(:));
IO.save(['mc_grid_' tag], GRID);
IO.save(['Y_smi_' tag], SMI.get_even_SH(eval_dirs, LMAX_FIT, CS));

% ------------------------------------------------------------- SMI, once
options = struct();
options.b     = bvals;
options.dirs  = bvecs;
options.sigma = sigma*ones(GRID);
options.mask  = true(GRID);
options.compartments  = {'IAS','EAS','FW'};
options.NoiseBias     = 'Rician';
options.Lmax          = [0 LMAX_FIT LMAX_FIT LMAX_FIT];
options.CS_phase      = CS;
options.D_FW          = D_FW;
options.flag_fit_fODF = 1;
options.fODF_regularization = struct('flag_nonneg',1,'lambda_nonneg',10, ...
                                     'lambda_tikhonov',0.3);
t0 = tic;
out = SMI.fit(dwi, options);
fprintf('SMI.fit (constrained) %.1f s\n', toc(t0));

L6  = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';
sc6 = sqrt((2*L6+1)/(4*pi))';

plmC = reshape(out.plm, [NVOX numel(L6)-1]);
shC  = [ones(NVOX,1) plmC].*repmat(sc6, NVOX, 1);
shC(~isfinite(shC)) = 0;
IO.save(['sh_smiC_' tag], shC);

% the same kernel, deconvolved without the constraint
s0 = out.RotInvs.S0(:,:,:,1);
t0 = tic;
plmU = SMI.get_plm_from_S_and_kernel(dwi./s0, options.Lmax, out.kernel, ...
        options.mask, bvals, ones(1,Ndwi), zeros(1,Ndwi), bvecs, CS, D_FW, []);
fprintf('unconstrained deconvolution %.1f s\n', toc(t0));
plmU = reshape(plmU, [NVOX numel(L6)-1]);
shU  = [ones(NVOX,1) plmU].*repmat(sc6, NVOX, 1);
shU(~isfinite(shU)) = 0;
IO.save(['sh_smiU_' tag], shU);

IO.save(['kernel_' tag], reshape(out.kernel, [NVOX size(out.kernel,4)]));
IO.save(['reginfo_' tag], [out.fODF_regularization.Niterations(:), ...
                           out.fODF_regularization.Nnegative_dirs(:), ...
                           out.fODF_regularization.flag_converged(:)]);
fprintf('wrote sh_smiC_%s and sh_smiU_%s [%d x %d]\n', tag, tag, NVOX, size(shC,2));
end

% =====================================================================
function G = pick_grid(N)
% A 3D grid holding exactly N voxels with every dimension > 1: SMI.vectorize
% takes a different branch if any spatial dimension is a singleton
% (README for Claude, section 4).
d = divisors_of(N);
d = d(d > 1 & d < N);
best = [];
for a = d
    if mod(N,a), continue, end
    m = N/a;
    e = divisors_of(m);
    e = e(e > 1 & e < m);
    for b = e
        c = m/b;
        if c > 1 && mod(m,b) == 0
            cand = sort([a b c]);
            if isempty(best) || (max(cand)-min(cand)) < (max(best)-min(best))
                best = cand;
            end
        end
    end
end
if isempty(best)
    error('pick_grid: cannot factor %d into three factors > 1', N);
end
G = best;
end

% =====================================================================
function d = divisors_of(n)
d = 1:floor(sqrt(n));
d = d(mod(n,d) == 0);
d = unique([d n./d]);
end

% =====================================================================
function m = rotate_about(n, deg)
% A unit vector at `deg` degrees from n, in an arbitrary but fixed plane.
n = n(:)'/norm(n);
t = [0 0 1]; if abs(n*t') > 0.9, t = [1 0 0]; end
e = t - (t*n')*n; e = e/norm(e);
m = cosd(deg)*n + sind(deg)*e;
m = m/norm(m);
end
