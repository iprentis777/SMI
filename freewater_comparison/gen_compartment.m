function gen_compartment(SNR, tag)
% gen_compartment(SNR, tag)
%
% Edema modelled as a REDISTRIBUTION of the non-intra-axonal water rather than
% as dilution. The intra-axonal fraction is held at f = 0.60 throughout and the
% remaining 0.40 is progressively converted from hindered extra-axonal water to
% free water:
%
%     (f, f_extra, fw) = (0.60, 0.40, 0.00)   healthy
%                        (0.60, 0.30, 0.10)
%                        (0.60, 0.20, 0.20)   <- the case asked for
%                        (0.60, 0.10, 0.30)
%                        (0.60, 0.00, 0.40)   all extracellular water freed
%
% This is a different perturbation from REPORT_fODF_freewater.md, where adding
% free water scaled f down with everything else. Here f never moves, so any
% weight built on f cannot dim these voxels at all -- which is exactly the
% thing worth measuring.
%
% Also writes two synthetic "brain" populations, healthy and edematous, from
% which CSD and MSMT-CSD can estimate their response functions the way the
% real tools do: from the most anisotropic voxels in the volume.

more off
IO = binio();

bvals     = IO.load('bvals')(:)';
bvecs     = IO.load('bvecs');
eval_dirs = IO.load('eval_dirs');
Ndwi      = numel(bvals);

LMAX_FIT = 6;
LMAX_GT  = 8;
CS       = 1;
D_FW     = 3;
NREP     = 40;
DA = 2.0; DEPAR = 2.0; DEPERP = 0.50; KAPPA = 16; F_INTRA = 0.60;

H  = fODF_modulation_helpers();
dq = H.dirs(2000);

n1 = [0.30 -0.50 0.81]; n1 = n1/norm(n1);
n2 = rotate_about(n1, 60);

FW_LIST = [0.00 0.10 0.20 0.30 0.40];
C = {};
for ig = 1:2
    if ig == 1, gname = 'single'; axes_ = {n1}; else, gname = 'cross60'; axes_ = {n1, n2}; end
    fod = zeros(size(dq,1),1);
    for k = 1:numel(axes_), fod = fod + H.watson(dq, axes_{k}, KAPPA); end
    for fw = FW_LIST
        s = struct();
        s.name   = sprintf('%s_fe%02d', gname, round((1-F_INTRA-fw)*100));
        s.fodf   = fod;
        s.axes   = axes_;
        s.kernel = [F_INTRA DA DEPAR DEPERP fw];   % f_extra = 1 - f - fw
        C{end+1} = s;
    end
end
s = struct();
s.name = 'csf'; s.fodf = ones(size(dq,1),1); s.axes = {};
s.kernel = [0.02 2.0 3.0 3.00 0.95];
C{end+1} = s;

NCOND = numel(C);
NVOX  = NCOND*NREP;
GRID  = [10 11 4];
assert(prod(GRID)==NVOX);

fprintf('compartment set: %d conditions x %d reps = %d voxels\n', NCOND, NREP, NVOX);
[S_clean, gt_amp, gt_amp6, cond_id, kern_gt, ax] = ...
    build_set(C, NREP, H, dq, bvals, bvecs, eval_dirs, LMAX_GT, LMAX_FIT, CS, D_FW);

sigma = 1/SNR;
rand('seed',777); randn('seed',777);
S_noisy = sqrt((S_clean + sigma*randn(size(S_clean))).^2 + ...
               (        sigma*randn(size(S_clean))).^2);

IO.save(['dwi_' tag], reshape(S_noisy, [GRID Ndwi]));
IO.save(['S_clean_' tag], S_clean);
IO.save('cond_id_c', cond_id);
IO.save('kern_gt_c', kern_gt);
IO.save('gt_amp_c', gt_amp);
IO.save('gt_amp6_c', gt_amp6);
IO.save('gt_axes_c', ax);
fid = fopen(fullfile(IO.dir(),'cond_names_c.txt'),'w');
for ic = 1:NCOND, fprintf(fid,'%s\n',C{ic}.name); end
fclose(fid);

% ------------------------------------------------------------- SMI fit
fprintf('fitting SMI (regularized, no modulation), SNR %g\n', SNR);
options = struct();
options.b = bvals; options.dirs = bvecs;
options.sigma = sigma*ones(GRID); options.mask = true(GRID);
options.compartments = {'IAS','EAS','FW'};
options.NoiseBias = 'Rician';
options.Lmax = [0 LMAX_FIT LMAX_FIT LMAX_FIT];
options.CS_phase = CS; options.D_FW = D_FW;
options.flag_fit_fODF = 1;
options.fODF_regularization = struct('flag_nonneg',1,'lambda_nonneg',10, ...
                                     'lambda_tikhonov',0.3);
t0 = tic;
out = SMI.fit(reshape(S_noisy,[GRID Ndwi]), options);
fprintf('SMI.fit took %.1f s\n', toc(t0));

plm = reshape(out.plm, [NVOX size(out.plm,4)]);
Y   = SMI.get_even_SH(eval_dirs, LMAX_FIT, CS);
L   = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';
flm = [ones(NVOX,1) plm] .* repmat(sqrt((2*L+1)/(4*pi))', NVOX, 1);
flm(~isfinite(flm)) = 0;
IO.save(['smi_amp_' tag], Y*flm');
IO.save(['smi_kernel_' tag], reshape(out.kernel, [NVOX size(out.kernel,4)]));
IO.save(['smi_pl_' tag], reshape(out.pl, [NVOX size(out.pl,4)]));
fprintf('wrote smi_amp_%s\n', tag);
end

% =====================================================================
function [S_clean, gt_amp, gt_amp6, cond_id, kern_gt, ax] = ...
    build_set(C, NREP, H, dq, bvals, bvecs, eval_dirs, LMAX_GT, LMAX_FIT, CS, D_FW)
NCOND = numel(C); NVOX = NCOND*NREP; Ndwi = numel(bvals);
S_clean = zeros(NVOX, Ndwi);
gt_amp  = zeros(size(eval_dirs,1), NCOND);
gt_amp6 = zeros(size(eval_dirs,1), NCOND);
cond_id = zeros(NVOX,1); kern_gt = zeros(NVOX,5);
Y_gt = SMI.get_even_SH(eval_dirs, LMAX_GT, CS);
L_gt = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';
keep6 = L_gt <= LMAX_FIT;
for ic = 1:NCOND
    plm_gt = H.mixture_plm(C{ic}.fodf, dq, LMAX_GT, CS);
    kv = [C{ic}.kernel 1 1];
    s  = H.signal(plm_gt, kv, bvals, ones(1,Ndwi), zeros(1,Ndwi), bvecs, ...
                  LMAX_GT, CS, D_FW);
    coef = [1; plm_gt(:)].*sqrt((2*L_gt+1)/(4*pi));
    gt_amp(:,ic)  = Y_gt*coef;
    gt_amp6(:,ic) = Y_gt(:,keep6)*coef(keep6);
    for ir = 1:NREP
        iv = (ic-1)*NREP + ir;
        S_clean(iv,:) = s(:)'; cond_id(iv) = ic; kern_gt(iv,:) = C{ic}.kernel;
    end
end
ax = nan(3,3,NCOND);
for ic = 1:NCOND
    for k = 1:numel(C{ic}.axes), ax(k,:,ic) = C{ic}.axes{k}; end
end
end

% =====================================================================
function m = rotate_about(n, deg)
n = n(:)'/norm(n);
t = [0 0 1]; if abs(n*t') > 0.9, t = [1 0 0]; end
e = t - (t*n')*n; e = e/norm(e);
m = cosd(deg)*n + sind(deg)*e; m = m/norm(m);
end
