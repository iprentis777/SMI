%% EXPERIMENTAL: what happens if p_00 is not fixed at 1?
% SMI's deconvolution imposes |p_00 = 1|, so its fODF integrates to 1 in every
% voxel and its isotropic floor is the constant |1/(4*pi) = 0.0796|. Two
% consequences are recorded in the handoff and both matter for tractography:
%
% * *the fODF carries no density information at all* -- a CSF voxel and a
%   coherent white matter voxel have equal mass;
% * *0.0796 sits above MRtrix's default |iFOD2 -cutoff 0.05|*, so an
%   unmodulated SMI fODF passes the termination test in EVERY voxel of the
%   brain, CSF included. Section 6.1 measured SMI leaving CSF at 0.32 where
%   MSMT-CSD leaves it at 0.028.
%
% This file asks what happens if that constraint is lifted and every
% coefficient including l = 0 is estimated from the data.
%
% *The comparison is one solver run two ways.* |helpers/fODF_free_l0_deconv.m|
% implements both conventions in one function, so the only thing that differs
% between the two arms is the convention itself. Step 2 checks that its
% |'fixed'| mode reproduces |SMI.get_plm_from_S_and_kernel| to machine
% precision -- if that check fails, nothing below is worth reading, because the
% |'free'| numbers would be coming from an unvalidated reimplementation.
%
% *The regularizer had to change, and in one substantive way.* With |p_00 = 1|
% the non-negativity constraint is INHOMOGENEOUS: the fixed floor |1/(4*pi)|
% moves to the right-hand side, and the fODF cannot be pushed below it. With
% |p_00| free the constraint is HOMOGENEOUS -- there is no floor, and the
% solution is free to shrink toward zero. That is not a detail of the
% implementation, it is the mechanism by which a free |p_00| could let SMI dim.
% Two smaller changes follow: the amplitude threshold |tau*mean(fODF)| is no
% longer a constant and is recomputed from the current solution, and Tikhonov
% damping never touches |l = 0| since shrinking the density term would defeat
% the purpose.
%
% *Three conditions, because the interesting one is not the easy one.* On
% healthy white matter deconvolved with its own kernel, |p_00| should come back
% at 1 and the two arms should agree -- that is the control, not the result.
% The question is what happens where the kernel does NOT describe the tissue.
%
%   1  healthy WM, 60 deg crossing, healthy kernel   matched: the control
%   2  edema signal, 60 deg crossing, healthy kernel mismatched, as on real data
%   3  pure free water, healthy kernel               the CSF termination case
%
% Condition 3 is the one to read first. It has no fibres, so there is no
% angular error to report -- the only question is whether the recovered fODF is
% dim enough for a tractography algorithm to stop.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
pkgdir = fileparts(here);
run(fullfile(pkgdir, 'oct_path.m'));
if exist('OCTAVE_VERSION', 'builtin'), warning('off', 'all'); end

H  = fODF_modulation_helpers();
RH = SMI_response_helpers();
FL = fODF_free_l0_deconv();
MC = mc_config();
VERDICT = {'** FAILED **', 'ok'};

%% Configuration
K_HEALTHY = [0.60 2.0 2.0 0.50 0.02];
K_EDEMA   = [0.10 2.4 2.7 1.15 0.35];
K_DECONV  = K_HEALTHY;      % the kernel handed to the deconvolution, always
D_FW      = 3;
KAPPA     = 16;
CROSS     = 60;
AXIS1     = [0.30 -0.50 0.81]; AXIS1 = AXIS1/norm(AXIS1);
LMAX_GT   = 8;
LMAX_FIT  = 6;
CS_PHASE  = 0;
PROTOCOL  = 'hcp_real_3shell.txt';
B0_SNAP   = 0.05;
NDIR_Q    = 3000;
NDIR_E    = 2000;
SEED      = 31415;
SEED_ORI  = 101;

SMOKE_TEST = true;
if SMOKE_TEST
    NORIENT = 4; NREP = 8;  SNR_LIST = [30 Inf];
else
    NORIENT = 8; NREP = 25; SNR_LIST = [10 30 Inf];
end

REG = struct('flag_nonneg',1, 'lambda_nonneg',1, 'tau',0.1, 'lambda_tikhonov',0);
IFOD2_CUTOFF = 0.05;        % MRtrix's default tractography termination threshold

fprintf('\n=== EXPERIMENTAL: free l = 0 deconvolution ===\n');
if SMOKE_TEST, fprintf('*** SMOKE_TEST = true: INDICATIVE ONLY ***\n'); end
fprintf('deconvolution kernel %s, Lmax %d, kappa %g\n\n', ...
        mat2str(K_DECONV), LMAX_FIT, KAPPA);

%% Step 1 -- protocol, geometry, signals for the three conditions
[bvals, bvecs] = MC.load_protocol_file(PROTOCOL);
bvals(bvals < B0_SNAP) = 0;
Ndwi = numel(bvals);

dq = H.dirs(NDIR_Q);
de = H.dirs(NDIR_E);
Lf = LMAX_FIT; nc = (Lf/2+1)*(Lf+1);
Lv = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
N_l = sqrt((2*Lv'+1)*(4*pi));
Y_dwi = SMI.get_even_SH(bvecs, Lf, CS_PHASE);
Ye    = SMI.get_even_SH(de,    Lf, CS_PHASE);
Yamp  = Ye .* repmat(sqrt((2*Lv'+1)/(4*pi)), size(Ye,1), 1);

% the design matrix, built exactly as SMI builds it
Kell = RH.Kell(K_DECONV, bvals, Lf, D_FW);                 % [Ndwi x Lf/2+1]
Kmat = Kell(:, repelem(1:(Lf/2+1), 2*(0:2:Lf)+1));         % [Ndwi x nc]
A    = Kmat .* (Y_dwi .* repmat(N_l, Ndwi, 1));

% non-negativity directions, as SMI picks them
dnn  = SMI.GetUniformHemisphereDirs(300);
Ynn  = SMI.get_even_SH(dnn, Lf, CS_PHASE);
Ynn  = Ynn .* repmat(sqrt((2*Lv'+1)/(4*pi)), size(Ynn,1), 1);

AX = cell(NORIENT,2);
AX{1,1} = AXIS1; AX{1,2} = MC.rotate_about(AXIS1, CROSS);
rand('state', SEED_ORI); randn('state', SEED_ORI);
for r = 2:NORIENT
    [Q,R0] = qr(randn(3)); Q = Q*diag(sign(diag(R0)));
    if det(Q) < 0, Q(:,1) = -Q(:,1); end
    AX{r,1} = (Q*AX{1,1}(:))'; AX{r,2} = (Q*AX{1,2}(:))';
end

COND = {struct('name','healthy WM, matched  ', 'kind','fib', 'K', K_HEALTHY), ...
        struct('name','edema, healthy kernel', 'kind','fib', 'K', K_EDEMA), ...
        struct('name','free water (CSF)     ', 'kind','iso', 'K', [])};

SIG = zeros(numel(COND), NORIENT, Ndwi);
for ic = 1:numel(COND)
    for r = 1:NORIENT
        if strcmp(COND{ic}.kind,'iso')
            SIG(ic,r,:) = exp(-bvals(:)*D_FW);        % pure isotropic, no fibres
        else
            f = H.watson(dq, AX{r,1}, KAPPA) + H.watson(dq, AX{r,2}, KAPPA);
            p = H.mixture_plm(f, dq, LMAX_GT, CS_PHASE);
            s = H.signal(p(:)', [COND{ic}.K 1 1], bvals, ones(1,Ndwi), ...
                         zeros(1,Ndwi), bvecs, LMAX_GT, CS_PHASE, D_FW);
            SIG(ic,r,:) = s(:);
        end
    end
end
fprintf('Step 1: %d conditions, %d orientations, %d volumes\n', numel(COND), NORIENT, Ndwi);
fprintf('        mean signal at b = 3 : ');
hib = bvals > 2.5;
for ic = 1:numel(COND)
    fprintf('%s %.4f   ', strtrim(COND{ic}.name), mean(squeeze(SIG(ic,1,hib))));
end
fprintf('\n\n');

%% Step 2 -- validate 'fixed' mode against SMI itself
% The whole experiment rests on this. If the reimplementation of the FIXED
% convention does not reproduce SMI's own output, the FREE numbers are
% meaningless.
regF = FL.reg_defaults(Lf, REG);
G3 = [2 2 2]; NV3 = prod(G3);
k4 = zeros([G3 5]); for j = 1:5, k4(:,:,:,j) = K_DECONV(j); end
sref = squeeze(SIG(1,1,:))';
plm_smi = SMI.get_plm_from_S_and_kernel(reshape(repmat(sref,NV3,1),[G3 Ndwi]), ...
            [0 Lf Lf Lf], k4, true(G3), bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
            bvecs, CS_PHASE, D_FW, REG);
p_ref = [1; squeeze(plm_smi(1,1,1,:))];
p_mine = FL.solve(A, sref(:), Ynn, regF, 'fixed');
e_val = max(abs(p_ref - p_mine));
fprintf('Step 2: CHECK ''fixed'' mode reproduces SMI.get_plm_from_S_and_kernel\n');
fprintf('        max|err| = %.2e   %s\n', e_val, VERDICT{1+(e_val < 1e-9)});
if e_val >= 1e-9
    error(['The fixed-convention solver does not match SMI. The free-l0 ' ...
           'results below would be from an unvalidated reimplementation.']);
end
fprintf('        (so a difference below is the CONVENTION, not the solver)\n\n');

%% Step 3 -- run both conventions on all three conditions
NSNR = numel(SNR_LIST); SGL = 1./SNR_LIST;
SNRLAB = cell(1,NSNR);
for is = 1:NSNR
    if isinf(SNR_LIST(is)), SNRLAB{is}='inf'; else, SNRLAB{is}=sprintf('%g',SNR_LIST(is)); end
end
MODES = {'fixed','free'};
P00 = nan(numel(COND),NSNR,2,NORIENT*NREP);
PKA = nan(numel(COND),NSNR,2,NORIENT*NREP);   % peak fODF amplitude (absolute)
AER = nan(numel(COND),NSNR,2,NORIENT*NREP);   % angular error, fibre conditions

fprintf('Step 3: deconvolving\n');
for ic = 1:numel(COND)
  for is = 1:NSNR
    rand('state', SEED+is); randn('state', SEED+is);
    for im = 1:2
      cnt = 0;
      for r = 1:NORIENT
        s0 = squeeze(SIG(ic,r,:))';
        for rep = 1:NREP
          cnt = cnt + 1;
          sg = SGL(is);
          if sg > 0
            sn = sqrt((s0 + sg*randn(1,Ndwi)).^2 + (sg*randn(1,Ndwi)).^2);
          else
            sn = s0;
          end
          p = FL.solve(A, sn(:), Ynn, regF, MODES{im});
          if any(~isfinite(p)), continue, end
          P00(ic,is,im,cnt) = p(1);
          amp = Yamp*p;
          [mx, imx] = max(amp);
          PKA(ic,is,im,cnt) = mx;
          if strcmp(COND{ic}.kind,'fib')
            v = de(imx,:);
            AER(ic,is,im,cnt) = min([acosd(min(abs(v*AX{r,1}(:)),1)), ...
                                     acosd(min(abs(v*AX{r,2}(:)),1))]);
          end
        end
      end
    end
  end
  fprintf('   %s done\n', strtrim(COND{ic}.name));
end

%% Step 4 -- results
fprintf('\n================ RESULTS ================\n');
fprintf('p00   the recovered l = 0 coefficient. FIXED mode pins it at 1 by\n');
fprintf('      construction; in FREE mode it is estimated and is the density.\n');
fprintf('peak  maximum fODF amplitude. MRtrix iFOD2 terminates below %.2f.\n', IFOD2_CUTOFF);
fprintf('pass  %% of voxels whose peak exceeds that cutoff, i.e. would be\n');
fprintf('      tracked THROUGH rather than terminated in.\n');
fprintf('err   angular error of the primary peak (fibre conditions only).\n\n');

for ic = 1:numel(COND)
    fprintf('--- %s ---\n', strtrim(COND{ic}.name));
    fprintf('  %-6s %-6s %9s %9s %9s %9s\n','SNR','mode','p00','peak','pass%','err');
    for is = 1:NSNR
        for im = 1:2
            p0 = squeeze(P00(ic,is,im,:)); p0 = p0(isfinite(p0));
            pk = squeeze(PKA(ic,is,im,:)); pk = pk(isfinite(pk));
            ae = squeeze(AER(ic,is,im,:)); ae = ae(isfinite(ae));
            if isempty(pk), continue, end
            if isempty(ae), es = '     --'; else, es = sprintf('%7.2f', mean(ae)); end
            fprintf('  %-6s %-6s %9.4f %9.4f %8.1f%% %s\n', SNRLAB{is}, MODES{im}, ...
                    mean(p0), mean(pk), 100*mean(pk > IFOD2_CUTOFF), es);
        end
    end
    fprintf('\n');
end

fprintf('================ THE NUMBER THIS FILE EXISTS FOR ================\n');
fprintf('Free water, noise free: can a tractography algorithm stop there?\n');
ABV = {'below','ABOVE'};        % cell index, not a local function: see the header
ii = find(isinf(SNR_LIST),1);
if ~isempty(ii)
    for im = 1:2
        pk = squeeze(PKA(3,ii,im,:)); pk = pk(isfinite(pk));
        p0 = squeeze(P00(3,ii,im,:)); p0 = p0(isfinite(p0));
        fprintf('  %-6s  p00 %7.4f   peak %7.4f   %s the %.2f cutoff\n', ...
                MODES{im}, mean(p0), mean(pk), ...
                ABV{1+(mean(pk) > IFOD2_CUTOFF)}, IFOD2_CUTOFF);
    end
end
fprintf('\n=== experiment done ===\n');
