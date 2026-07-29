function gen_and_fit_smi(SNR, tag)
% gen_and_fit_smi(SNR, tag)
%
% Builds the ground truth, synthesises the noisy DWI, and fits it with the
% REGULARIZED SMI deconvolution (no modulation). Writes the DWI out so that
% dipy's CSD and MSMT-CSD see byte-identical data, and writes SMI's fODF
% sampled on the shared evaluation sphere so that peak extraction downstream
% is the same code for all three methods.
%
% Ground truth is a WM tissue kernel (stick + zeppelin) plus a free water
% compartment. Note SMI's convention, SMI.m:2075:  f_extra = 1 - f - fw,
% i.e. f, f_extra and fw are all fractions of the TOTAL signal. So holding
% the tissue microstructure fixed while adding free water means
%     f = f_intra_tissue * (1 - fw)
% which is what FW_LIST below does.

more off
IO = binio();

% ------------------------------------------------------------------ setup
bvals     = IO.load('bvals')(:)';           % [1 x Ndwi], ms/um^2
bvecs     = IO.load('bvecs');               % [Ndwi x 3]
eval_dirs = IO.load('eval_dirs');           % [Ndir x 3]
Ndwi      = numel(bvals);

LMAX_FIT = 6;                 % all three methods fit at this angular order
LMAX_GT  = 8;                 % ground truth signal is NOT band limited to the
                              % fit order. 8 is the ceiling: the SM kernel
                              % hard-codes Legendre polynomials only to l=8
                              % (SMI.m:2099-2107).
CS       = 1;                 % SMI.fit's default CS_phase
D_FW     = 3;                 % um^2/ms
NREP     = 40;

% WM tissue microstructure, held fixed across free water levels
F_INTRA = 0.60;  DA = 2.0;  DEPAR = 2.0;  DEPERP = 0.50;
KAPPA   = 16;

H = fODF_modulation_helpers();
dq = H.dirs(2000);                          % quadrature grid for the fODFs

% fibre geometry
n1 = [0.30 -0.50 0.81]; n1 = n1/norm(n1);
n2_60 = rotate_about(n1, 60);
n2_45 = rotate_about(n1, 45);

% ------------------------------------------------------------- conditions
% {name, fODF amplitudes on dq, true peak axes, kernel [f Da Depar Deperp fw]}
FW_LIST = [0.0 0.4];
C = {};
geoms = { 'single',  {n1}, ...
          'cross60', {n1, n2_60}, ...
          'cross45', {n1, n2_45} };
for ig = 1:3
    gname = geoms{2*ig-1};
    axes_ = geoms{2*ig};
    fod = zeros(size(dq,1),1);
    for k = 1:numel(axes_)
        fod = fod + H.watson(dq, axes_{k}, KAPPA);
    end
    for fw = FW_LIST
        % built field by field: struct('axes',cell) would make a struct ARRAY
        s = struct();
        s.name   = sprintf('%s_fw%02d', gname, round(fw*100));
        s.fodf   = fod;
        s.axes   = axes_;
        s.kernel = [F_INTRA*(1-fw) DA DEPAR DEPERP fw];
        C{end+1} = s;
    end
end
% pure CSF, for the scale-free WM/CSF contrast metric
s = struct();
s.name   = 'csf';
s.fodf   = ones(size(dq,1),1);
s.axes   = {};
s.kernel = [0.02 2.0 3.0 3.00 0.95];
C{end+1} = s;

NCOND = numel(C);
NVOX  = NCOND*NREP;
% SMI.vectorize takes a different branch if any spatial dim is a singleton
% (README section 4), so all three dimensions must be > 1.
GRID = [10 14 2];
assert(prod(GRID)==NVOX, 'grid %s does not hold %d voxels', mat2str(GRID), NVOX);

% ------------------------------------------------- ground truth + signals
fprintf('building ground truth (%d conditions x %d reps = %d voxels)\n', ...
        NCOND, NREP, NVOX);
S_clean  = zeros(NVOX, Ndwi);
gt_amp   = zeros(size(eval_dirs,1), NCOND);
gt_amp6  = zeros(size(eval_dirs,1), NCOND);
cond_id  = zeros(NVOX,1);
kern_gt  = zeros(NVOX,5);
Y_gt     = SMI.get_even_SH(eval_dirs, LMAX_GT, CS);
L_gt     = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';

for ic = 1:NCOND
    plm_gt = H.mixture_plm(C{ic}.fodf, dq, LMAX_GT, CS);
    kv     = [C{ic}.kernel 1 1];                       % [f Da Depar Deperp fw T2a T2e]
    s      = H.signal(plm_gt, kv, bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                      bvecs, LMAX_GT, CS, D_FW);
    % ground truth fODF amplitude on the shared evaluation sphere
    gt_amp(:,ic) = Y_gt*([1; plm_gt(:)].*sqrt((2*L_gt+1)/(4*pi)));
    % and the same truncated to the fit order: the ceiling any Lmax 6 method
    % could reach, which is the honest reference for angular resolution
    keep6 = L_gt <= LMAX_FIT;
    gt_amp6(:,ic) = Y_gt(:,keep6)*([1; plm_gt(:)](keep6).*sqrt((2*L_gt(keep6)+1)/(4*pi)));
    for ir = 1:NREP
        iv = (ic-1)*NREP + ir;
        S_clean(iv,:) = s(:)';
        cond_id(iv)   = ic;
        kern_gt(iv,:) = C{ic}.kernel;
    end
end

% Rician noise. S0 = 1 by construction, so sigma = 1/SNR.
sigma = 1/SNR;
rand('seed',12345); randn('seed',12345);
S_noisy = sqrt((S_clean + sigma*randn(size(S_clean))).^2 + ...
               (        sigma*randn(size(S_clean))).^2);

dwi = reshape(S_noisy, [GRID Ndwi]);
IO.save(['dwi_' tag], dwi);
IO.save(['S_clean_' tag], S_clean);
IO.save('cond_id', cond_id);
IO.save('kern_gt', kern_gt);
IO.save('gt_amp', gt_amp);
IO.save('gt_amp6', gt_amp6);
IO.save('grid', GRID);
% true fibre axes, padded to 3 rows per condition (NaN where absent)
ax = nan(3,3,NCOND);
for ic = 1:NCOND
    for k = 1:numel(C{ic}.axes), ax(k,:,ic) = C{ic}.axes{k}; end
end
IO.save('gt_axes', ax);
fid = fopen(fullfile(IO.dir(),'cond_names.txt'),'w');
for ic = 1:NCOND, fprintf(fid,'%s\n',C{ic}.name); end
fclose(fid);

% ------------------------------------------------------------- SMI fit
fprintf('fitting SMI (regularized deconvolution, no modulation), SNR %g\n', SNR);
options = struct();
options.b     = bvals;
options.dirs  = bvecs;
options.sigma = sigma*ones(GRID);
options.mask  = true(GRID);
options.compartments = {'IAS','EAS','FW'};   % default is IAS+EAS only
options.NoiseBias    = 'Rician';
options.Lmax         = [0 LMAX_FIT LMAX_FIT LMAX_FIT];
options.CS_phase     = CS;
options.D_FW         = D_FW;
options.flag_fit_fODF = 1;
options.fODF_regularization = struct('flag_nonneg',1,'lambda_nonneg',10, ...
                                     'lambda_tikhonov',0.3);
% modulation deliberately left off (flag_modulate defaults to 0)

t0 = tic;
out = SMI.fit(dwi, options);
fprintf('SMI.fit took %.1f s\n', toc(t0));

% ---------------------------------- SMI fODF on the shared eval sphere
plm = reshape(out.plm, [NVOX size(out.plm,4)]);
Y   = SMI.get_even_SH(eval_dirs, LMAX_FIT, CS);
L   = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';
flm = [ones(NVOX,1) plm] .* repmat(sqrt((2*L+1)/(4*pi))', NVOX, 1);
flm(~isfinite(flm)) = 0;
smi_amp = Y*flm';                                  % [Ndir x NVOX]

IO.save(['smi_amp_' tag], smi_amp);
IO.save(['smi_kernel_' tag], reshape(out.kernel, [NVOX size(out.kernel,4)]));
IO.save(['smi_pl_' tag], reshape(out.pl, [NVOX size(out.pl,4)]));
fprintf('wrote smi_amp_%s [%d x %d]\n', tag, size(smi_amp,1), size(smi_amp,2));
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
