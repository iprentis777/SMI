function sweep_nonneg(SNR, NREP, tag)
% sweep_nonneg(SNR, NREP, tag)
%
% The non-negativity weight, swept on the same crossing conditions as
% gen_montecarlo.m.
%
% This exists because the main comparison found that the SHIPPED default,
% lambda_nonneg = 10, does not resolve the 45 degree crossing at all, while the
% unconstrained deconvolution of the SAME kernel resolves it almost every time.
% That reproduces the observation recorded in "README for Claude" section 2
% item 3, which was the most concrete open task in the repository, and this
% sweep is what turns it from an observation into a number per weight.
%
% One SMI.fit, then one deconvolution per setting from the same kernel and the
% same normalized signal, so the arms differ only in the regularizer.
%
% Writes sh_sweep<k>_<tag> for k = 1..numel(SETTINGS) and a text file naming
% them, which score_sweep.py reads.

more off
IO = binio();
MR = mrtrix_io();
mdir = fullfile(fileparts(mfilename('fullpath')), 'mrtrix');
if ~exist(mdir,'dir'), mkdir(mdir); end

eval_dirs = IO.load('eval_dirs');

% Same constants as gen_montecarlo.m, from the same place, so the sweep and the
% main run cannot disagree about what is being simulated. Only the RNG seed
% differs (C.SEED_SWEEP), which keeps the two sets of realisations independent.
C        = mc_config();
[bvals, bvecs] = C.load_protocol();
Ndwi     = numel(bvals);
LMAX_FIT = C.LMAX_FIT;
LMAX_GT  = C.LMAX_GT;
CS       = C.CS_PHASE;
D_FW     = C.D_FW;
KAPPA    = C.KAPPA;
ANGLES   = C.ANGLES;
K_WM     = C.K_WM;

% {label, lambda_nonneg (0 = constraint off), lambda_tikhonov}
SETTINGS = { 'nonneg off, tik 0.3',   0, 0.3
             'nonneg 1,   tik 0.3',   1, 0.3
             'nonneg 3,   tik 0.3',   3, 0.3
             'nonneg 10,  tik 0.3',  10, 0.3
             'nonneg 10,  tik 0',    10, 0.0
             'nonneg 3,   tik 0',     3, 0.0 };

H  = fODF_modulation_helpers();
dq = H.dirs(C.NDIR_Q);

NCOND = numel(ANGLES);
NVOX  = NCOND*NREP;
GRID  = C.pick_grid(NVOX);
fprintf('sweep: SNR %g, %d x %d = %d voxels, grid %s\n', ...
        SNR, NCOND, NREP, NVOX, mat2str(GRID));

L_gt  = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';
keep6 = L_gt <= LMAX_FIT;
S_cond = zeros(NCOND, Ndwi);
sh_gt6 = zeros(NCOND, sum(keep6));
ax     = nan(2, 3, NCOND);
for ic = 1:NCOND
    axes_ = C.condition_axes(ic);
    fod = zeros(size(dq,1),1);
    for k = 1:numel(axes_), fod = fod + H.watson(dq, axes_{k}, KAPPA); end
    plm_gt = H.mixture_plm(fod, dq, LMAX_GT, CS);
    coef   = [1; plm_gt(:)].*sqrt((2*L_gt+1)/(4*pi));
    sh_gt6(ic,:) = coef(keep6)';
    s = H.signal(plm_gt, [K_WM 1 1], bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                 bvecs, LMAX_GT, CS, D_FW);
    S_cond(ic,:) = s(:)';
    for k = 1:numel(axes_), ax(k,:,ic) = axes_{k}; end
end

cond_id = repelem((1:NCOND)', NREP, 1);
sigma   = 1/SNR;
rand('seed', C.SEED_SWEEP); randn('seed', C.SEED_SWEEP);
S_clean = S_cond(cond_id, :);
S_noisy = sqrt((S_clean + sigma*randn(size(S_clean))).^2 + ...
               (         sigma*randn(size(S_clean))).^2);
dwi = reshape(S_noisy, [GRID Ndwi]);

IO.save(['mc_dwi_' tag], dwi);
IO.save(['mc_cond_id_' tag], cond_id);
IO.save(['mc_gt_axes_' tag], ax);
IO.save(['mc_sh_gt6_' tag], sh_gt6);
IO.save(['mc_angles_' tag], ANGLES(:));
IO.save(['Y_smi_' tag], SMI.get_even_SH(eval_dirs, LMAX_FIT, CS));
MR.write(fullfile(mdir,['mask_' tag]), ones(GRID), struct('datatype','UInt8'));
gtv = zeros([GRID sum(keep6)]);
gtc = reshape(sh_gt6(cond_id,:), [GRID sum(keep6)]);
gtv(:) = gtc(:);
MR.write(fullfile(mdir,['gtfod_' tag]), gtv);

options = struct();
options.b = bvals; options.dirs = bvecs;
options.sigma = sigma*ones(GRID); options.mask = true(GRID);
options.compartments  = {'IAS','EAS','FW'};
options.NoiseBias     = 'Rician';
options.Lmax          = [0 LMAX_FIT LMAX_FIT LMAX_FIT];
options.CS_phase      = CS; options.D_FW = D_FW;
options.flag_fit_fODF = 0;                 % the kernel only; the fODFs below
t0 = tic;
out = SMI.fit(dwi, options);
fprintf('SMI.fit (kernel only) %.1f s\n', toc(t0));

s0  = out.RotInvs.S0(:,:,:,1);
dn  = dwi./s0;
L6  = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';
sc6 = sqrt((2*L6+1)/(4*pi))';

fid = fopen(fullfile(IO.dir(),['sweep_names_' tag '.txt']),'w');
for is = 1:size(SETTINGS,1)
    if SETTINGS{is,2} == 0
        reg = struct('flag_nonneg',0,'lambda_tikhonov',SETTINGS{is,3});
    else
        reg = struct('flag_nonneg',1,'lambda_nonneg',SETTINGS{is,2}, ...
                     'lambda_tikhonov',SETTINGS{is,3});
    end
    t0 = tic;
    plm = SMI.get_plm_from_S_and_kernel(dn, options.Lmax, out.kernel, ...
            options.mask, bvals, ones(1,Ndwi), zeros(1,Ndwi), bvecs, CS, D_FW, reg);
    plm = reshape(plm, [NVOX numel(L6)-1]);
    sh  = [ones(NVOX,1) plm].*repmat(sc6, NVOX, 1);
    sh(~isfinite(sh)) = 0;
    IO.save(sprintf('sh_sweep%d_%s', is, tag), sh);
    MR.write(fullfile(mdir,sprintf('sweep%dfod_%s', is, tag)), ...
             reshape(sh, [GRID numel(L6)]));
    fprintf('  %-22s %5.1f s\n', SETTINGS{is,1}, toc(t0));
    fprintf(fid, '%s\n', SETTINGS{is,1});
end
fclose(fid);
end

% pick_grid, divisors_of and rotate_about used to live here as local functions,
% duplicating gen_montecarlo.m's copies. They are mc_config.m's now.
