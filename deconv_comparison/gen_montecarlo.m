function gen_montecarlo(SNR, NREP, tag)
% gen_montecarlo(SNR, NREP, tag)
%
% The Monte Carlo experiment of Reports/REPORT_SMI_deconvolution_MonteCarlo.md, in the
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
% This file produces the SMI arm and the inputs the MRtrix arms need. CSD and
% MSMT-CSD are run by MRtrix3 itself (run_mrtrix.sh), not by a reimplementation,
% so everything they touch is written out in MRtrix's own image format.
%
% CS_phase = 0 THROUGHOUT. With SMI's default CS_phase = 1 the spherical
% harmonic basis differs from MRtrix's by (-1)^m, which is a 180 degree
% rotation about z of every fODF -- verified against MRtrix in
% check_mrtrix_basis.sh. Setting it to 0 makes SMI's coefficients literally
% MRtrix's, so no conversion sits between the two sides of the comparison.
%
% Written per SNR, in data/:
%   mc_dwi_<tag>            the noisy signal (binio, for reference)
%   sh_smi_<tag>            SMI fODF, SH coefficients
%   kernel_<tag>, s0_<tag>  the fitted kernel and the b=0 invariant
% and in mrtrix/:
%   mc_<tag>.mih            the same signal as an MRtrix image with dw_scheme
%   smifod_<tag>.mih        the SMI fODF, directly in MRtrix's SH basis
%   gtfod_<tag>.mih         the band limited ground truth, same basis
%   mask_<tag>.mih          all ones (per tag: the arms have different sizes)

more off
IO = binio();
MR = mrtrix_io();

bvals     = IO.load('bvals'); bvals = bvals(:)';
bvecs     = IO.load('bvecs');
Ndwi      = numel(bvals);

% Every constant of the experiment comes from mc_config.m, which sweep_nonneg.m
% reads too, so the main run and the regularizer sweep cannot disagree about
% what is being simulated.
C        = mc_config();
LMAX_FIT = C.LMAX_FIT;
LMAX_GT  = C.LMAX_GT;
CS       = C.CS_PHASE;
D_FW     = C.D_FW;
KAPPA    = C.KAPPA;
ANGLES   = C.ANGLES;
K_WM     = C.K_WM;

H  = fODF_modulation_helpers();
dq = H.dirs(C.NDIR_Q);

NCOND = numel(ANGLES);
NVOX  = NCOND*NREP;
GRID  = C.pick_grid(NVOX);
mdir  = fullfile(fileparts(mfilename('fullpath')), 'mrtrix');
if ~exist(mdir,'dir'), mkdir(mdir); end
fprintf('SNR %g, %d conditions x %d reps = %d voxels, grid %s\n', ...
        SNR, NCOND, NREP, NVOX, mat2str(GRID));

% --------------------------------------------------- ground truth signals
Y_gt  = SMI.get_even_SH(dq, LMAX_GT, CS);
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
S_clean = S_cond(cond_id, :);

% ------------------------------------------------------------ Rician noise
% S0 = 1 by construction, so sigma = 1/SNR. Complex Gaussian noise on a real
% signal, magnitude taken: exactly Rician.
sigma = 1/SNR;
rand('seed', C.SEED_MC); randn('seed', C.SEED_MC);
S_noisy = sqrt((S_clean + sigma*randn(size(S_clean))).^2 + ...
               (         sigma*randn(size(S_clean))).^2);
dwi = reshape(S_noisy, [GRID Ndwi]);

IO.save(['mc_dwi_' tag], dwi);
IO.save(['mc_cond_id_' tag], cond_id);
IO.save(['mc_gt_axes_' tag], ax);
IO.save(['mc_sh_gt6_' tag], sh_gt6);
IO.save(['mc_angles_' tag], ANGLES(:));
IO.save(['mc_grid_' tag], GRID);

% ------------------------------------------------------- MRtrix inputs
% b in s/mm^2, which is what MRtrix expects; the protocol is carried in
% ms/um^2 everywhere else in this package.
MR.write(fullfile(mdir,['mc_' tag]), dwi, struct('grad',[bvecs bvals(:)*1000]));
MR.write(fullfile(mdir,['mask_' tag]), ones(GRID), struct('datatype','UInt8'));
gt_vol = zeros([GRID sum(keep6)]);
gtc = reshape(sh_gt6(cond_id,:), [GRID sum(keep6)]);
gt_vol(:) = gtc(:);
MR.write(fullfile(mdir,['gtfod_' tag]), gt_vol);

% ------------------------------------------------------------- SMI
% The constrained deconvolution at the shipped defaults. lambda_nonneg is 1
% (SMI.m:978); the sweep behind that choice is section 6 of the report.
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
options.fODF_regularization = struct('flag_nonneg',1,'lambda_tikhonov',0.3);
t0 = tic;
out = SMI.fit(dwi, options);
fprintf('SMI.fit %.1f s (lambda_nonneg = %g)\n', toc(t0), ...
        out.fODF_regularization.lambda_nonneg);

L6  = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';
sc6 = sqrt((2*L6+1)/(4*pi))';
plm = reshape(out.plm, [NVOX numel(L6)-1]);
sh  = [ones(NVOX,1) plm].*repmat(sc6, NVOX, 1);
sh(~isfinite(sh)) = 0;

IO.save(['sh_smi_' tag], sh);
IO.save(['kernel_' tag], reshape(out.kernel, [NVOX size(out.kernel,4)]));
IO.save(['s0_' tag], out.RotInvs.S0(:,:,:,1));
IO.save(['reginfo_' tag], [out.fODF_regularization.Niterations(:), ...
                           out.fODF_regularization.Nnegative_dirs(:), ...
                           out.fODF_regularization.flag_converged(:)]);
MR.write(fullfile(mdir,['smifod_' tag]), reshape(sh, [GRID numel(L6)]));
fprintf('wrote sh_smi_%s [%d x %d] and mrtrix/smifod_%s.mih\n', ...
        tag, NVOX, size(sh,2), tag);
end

% pick_grid, divisors_of and rotate_about used to live here as local functions.
% They are mc_config.m's now -- keeping a second copy is how the walkthrough and
% the pipeline would have drifted apart.
