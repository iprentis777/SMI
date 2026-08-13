%% Parameter sweep: which deconvolution settings give the most accurate fODF?
% A formal sampling of the tunable parameter space of all three deconvolution
% arms, on one shared simulation, scored two different ways.
%
% *The kernel is ground truth AND is handed to every arm.* No kernel is
% estimated anywhere in this file -- SMI goes through
% |SMI.get_plm_from_S_and_kernel| with the true kernel, exactly as Arm 1 of
% |notebooks/smi_wm_60deg.m| does, and the CSD arms get the response derived
% from that same kernel by |r_l = K_l sqrt((2l+1)4pi)|. So a difference between
% two rows of the output is the SETTING and nothing else.
%
% *Two families of metric, reported separately and never averaged together.*
% This is the point of the file rather than a detail. The repository already
% contains two sweeps that disagree about |lambda_nonneg|:
% |examples/example_fODF_regularization_sweep.m| minimises relative L2 error
% over the sphere and prefers 10, while the Monte Carlo in |Reports/| scores
% peaks and prefers 1. They are not in conflict -- they score different things,
% and which one matters depends on what consumes the fODF.
%
%   PEAK metrics      what tractography consumes
%     err   mean angular error of the primary peak, from sh2peaks
%     cor   fraction of realisations recovering exactly the true fibre count
%     spur  spurious peaks per voxel
%
%   SHAPE metrics     what an fODF-domain analysis consumes
%     ACC   angular correlation over l >= 2 against the band-limited truth
%           at the SAME Lmax, so the band limit is not charged as method error
%     L2    relative L2 error over the sphere, both sides on p_00 = 1
%     neg   fraction of the sphere where the recovered fODF is negative
%
% Step 4 prints where the two families pick different winners. That
% disagreement is a result, not a defect in the sweep.
%
% *No local functions.* MATLAB requires a script's local functions at the END
% of the file while Octave cannot call them there at all, so every arm's fODFs
% are collected first and scored in one pass afterwards.
%
% *Runtime.* No |SMI.fit| anywhere -- the SMI arm is the fixed-kernel solve, so
% there is no regression training, and the MRtrix arms cost seconds each.
% |SMOKE_TEST| cuts the grids to a handful of points.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
pkgdir = fileparts(here);                    % deconv_comparison/
run(fullfile(pkgdir, 'oct_path.m'));
if exist('OCTAVE_VERSION', 'builtin'), warning('off', 'all'); end

H  = fODF_modulation_helpers();
RH = SMI_response_helpers();
MC = mc_config();
MR = mrtrix_io();
VERDICT = {'** FAILED **', 'ok'};

%% Configuration
% The simulated tissue and geometry are IDENTICAL to
% notebooks/smi_wm_60deg.m. Nothing here differs from that file except the
% settings being swept.

K_WM     = [0.60 2.0 2.0 0.50 0.02];  % ground truth AND handed to every arm
D_FW     = 3;
D_GM     = 0.8;                       % only to build msmt_csd's second tissue
KAPPA    = 16;
CROSS    = 60;
AXIS1    = [0.30 -0.50 0.81]; AXIS1 = AXIS1/norm(AXIS1);
LMAX_GT  = 8;
CS_PHASE = 0;
PROTOCOL = 'hcp_real_3shell.txt';
B0_SNAP  = 0.05;
NDIR_Q   = 3000;
SEED     = 31415;
SEED_ORI = 101;

SMOKE_TEST = true;
if SMOKE_TEST
    NORIENT = 4;  NREP = 10;
    SNR_LIST = [10 Inf];
    LMAX_LIST = 6;
else
    NORIENT = 12; NREP = 40;
    SNR_LIST = [10 30 Inf];
    LMAX_LIST = [4 6 8];
end

PEAK_REL = 0.30;
PEAK_NUM = 3;
NDIR_E   = 4000;    % directions for the SHAPE metrics

% Where the transient MRtrix images go. A full-size run writes a lot, and
% writing it into a synced folder is how a truncated image happens.
SCRATCH_DIR = '';
KEEP_FILES  = false;

%% The sweeps
% SMI knobs are the fields of SMI.fODF_RegularizationDefaults; MRtrix knobs are
% the documented options of dwi2fod. One-at-a-time around each default, plus
% one 2-D grid per arm on the pair that interacts.

SMI_BASE = struct('flag_nonneg',1, 'lambda_nonneg',1, 'tau',0.1, ...
                  'Niter',50, 'Ndirs',300, 'Lmax_init',4, ...
                  'max_neg_fraction',0.9, 'lambda_tikhonov',0);

if SMOKE_TEST
    SMI_OAT  = { {'lambda_nonneg',[0.3 3]}, {'lambda_tikhonov',[0.3]} };
    SMI_GRID = {'lambda_nonneg',[0.3 1 3], 'tau',[0.05 0.1]};
    CSD_NEG  = [1e-10 1];   CSD_NORM = [1e-10 1];   CSD_THR = 0;
    CSD_FILT = {''};
    MSMT_NEG = [1e-10 1];   MSMT_NORM = [1e-10 1e-3];
else
    SMI_OAT  = { {'lambda_nonneg',   [0.1 0.3 3 10]}, ...
                 {'tau',             [0.02 0.05 0.2 0.4]}, ...
                 {'lambda_tikhonov', [0.03 0.1 0.3 1]}, ...
                 {'Ndirs',           [100 200 600]}, ...
                 {'Lmax_init',       [2 6]}, ...
                 {'Niter',           [10 25 100]}, ...
                 {'max_neg_fraction',[0.5 0.7]} };
    SMI_GRID = {'lambda_nonneg',[0.1 0.3 1 3 10], 'tau',[0.02 0.05 0.1 0.2]};
    CSD_NEG  = [1e-10 1e-3 0.1 1 3];  CSD_NORM = [1e-10 1e-3 0.1 1];
    CSD_THR  = [0 0.05];              CSD_FILT = {'', '1,1,1,1,1'};
    MSMT_NEG = [1e-10 1e-3 0.1 1 3];  MSMT_NORM = [1e-10 1e-3 0.1 1];
end

% flag_nonneg = 0 is the unconstrained reference: what the shipped toolbox does
% by default, and the only rows here that are not constrained deconvolutions.
SMI_EXTRA = { struct('flag_nonneg',0, 'lambda_tikhonov',0), ...
              struct('flag_nonneg',0, 'lambda_tikhonov',0.3) };

fprintf('\n=== deconvolution parameter sweep ===\n');
if SMOKE_TEST, fprintf('*** SMOKE_TEST = true: reduced grids, INDICATIVE ONLY ***\n'); end
fprintf('kernel %s (ground truth AND handed to every arm), kappa %g, %g deg\n', ...
        mat2str(K_WM), KAPPA, CROSS);
fprintf('%d orientations x %d reps x %d SNR, Lmax %s\n\n', ...
        NORIENT, NREP, numel(SNR_LIST), mat2str(LMAX_LIST));

%% Step 1 -- protocol, ground truth, signal
[bvals, bvecs] = MC.load_protocol_file(PROTOCOL);
bvals(bvals < B0_SNAP) = 0;
Ndwi = numel(bvals);
[tbl, ~, shell_id] = SMI.Group_dwi_in_shells_b_beta_TE(bvals, [], [], []);
b_shell = tbl(1,:); n_shell = tbl(3,:); nsh = numel(b_shell);

dq   = H.dirs(NDIR_Q);
de   = H.dirs(NDIR_E);
L_gt = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';

AX = cell(NORIENT,2);
AX{1,1} = AXIS1; AX{1,2} = MC.rotate_about(AXIS1, CROSS);
rand('state', SEED_ORI); randn('state', SEED_ORI);
for r = 2:NORIENT
    [Q,R0] = qr(randn(3)); Q = Q*diag(sign(diag(R0)));
    if det(Q) < 0, Q(:,1) = -Q(:,1); end
    AX{r,1} = (Q*AX{1,1}(:))'; AX{r,2} = (Q*AX{1,2}(:))';
end
e_sep   = 0;
sh_gt   = zeros(NORIENT, numel(L_gt));
S_clean = zeros(NORIENT, Ndwi);
for r = 1:NORIENT
    e_sep = max(e_sep, abs(acosd(abs(AX{r,1}*AX{r,2}')) - CROSS));
    f = H.watson(dq, AX{r,1}, KAPPA) + H.watson(dq, AX{r,2}, KAPPA);
    p = H.mixture_plm(f, dq, LMAX_GT, CS_PHASE);
    sh_gt(r,:) = ([1; p(:)] .* sqrt((2*L_gt+1)/(4*pi)))';
    s = H.signal(p(:)', [K_WM 1 1], bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                 bvecs, LMAX_GT, CS_PHASE, D_FW);
    S_clean(r,:) = s(:)';
end
fprintf('Step 1: CHECK crossings are %g deg   max err %.2e   %s\n', ...
        CROSS, e_sep, VERDICT{1+(e_sep < 1e-9)});

NSNR = numel(SNR_LIST); SIG = 1./SNR_LIST;
SNRLAB = cell(1,NSNR);
for is = 1:NSNR
    if isinf(SNR_LIST(is)), SNRLAB{is} = 'inf';
    else, SNRLAB{is} = sprintf('%g', SNR_LIST(is)); end
end
NVOX_SNR = NORIENT*NREP; NVOX = NVOX_SNR*NSNR;
orient_id = repmat(repelem((1:NORIENT)', NREP, 1), NSNR, 1);
snr_id    = repelem((1:NSNR)', NVOX_SNR, 1);
S_rep   = S_clean(orient_id,:);
S_noisy = zeros(NVOX, Ndwi);
for is = 1:NSNR
    rows = find(snr_id == is); sg = SIG(is);
    rand('state', SEED+is); randn('state', SEED+is);
    Sb = S_rep(rows,:);
    S_noisy(rows,:) = sqrt((Sb + sg*randn(size(Sb))).^2 + (sg*randn(size(Sb))).^2);
end
GRID_SNR = MC.pick_grid(NVOX_SNR);
GRID_ALL = MC.pick_grid(NVOX);
fprintf('        %d voxels (%d per SNR block, grid %s)\n\n', NVOX, NVOX_SNR, mat2str(GRID_SNR));

if isempty(SCRATCH_DIR), sroot = pkgdir; else, sroot = SCRATCH_DIR; end
WDIR = fullfile(sroot, 'sweep_tmp');
if ~exist(WDIR, 'dir'), mkdir(WDIR); end

f_dwi = fullfile(WDIR,'dwi'); f_mask = fullfile(WDIR,'mask'); f_b3 = fullfile(WDIR,'b3');
MR.write([f_dwi '.mif'], reshape(S_noisy,[GRID_ALL Ndwi]), ...
         struct('grad',[bvecs bvals(:)*1000]));
MR.write([f_mask '.mif'], ones(GRID_ALL), struct('datatype','UInt8'));
[~,t1] = system(sprintf('mrinfo -shell_bvalues "%s.mif" 2>/dev/null', f_dwi));
b_mr = sscanf(t1,'%f')';
[st,tx] = system(sprintf('dwiextract "%s.mif" -shells 0,%g "%s.mif" -force -quiet 2>&1', ...
                         f_dwi, b_mr(end), f_b3));
if st ~= 0, fprintf(2,'%s\n',tx); error('dwiextract failed'); end
fprintf('        MRtrix shells %s s/mm^2\n\n', mat2str(round(b_mr)));

%% Step 2 -- run every setting, collect the fODFs
% Nothing is scored here. Each entry of JOB carries the arm, the setting name,
% the Lmax and the [NVOX x ncoef] SH matrix, and Step 3 scores them all through
% one identical path so a difference between rows cannot be a scoring
% difference.

JOB = {};
fprintf('Step 2: running the sweep\n');
for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);
    nc = (Lf/2+1)*(Lf+1);
    Lv = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
    sc = sqrt((2*Lv+1)/(4*pi));

    r_wm  = RH.zh(K_WM, b_mr/1000, Lf, D_FW);
    r_gm  = exp(-(b_mr(:)/1000)*D_GM)*sqrt(4*pi);
    r_csf = exp(-(b_mr(:)/1000)*D_FW)*sqrt(4*pi);
    Rw = fullfile(WDIR,sprintf('wm%d.txt',Lf));  RH.write_response(Rw, r_wm);
    Rg = fullfile(WDIR,'gm.txt');                RH.write_response(Rg, r_gm);
    Rc = fullfile(WDIR,'csf.txt');               RH.write_response(Rc, r_csf);
    Rb = fullfile(WDIR,sprintf('wb%d.txt',Lf));  RH.write_response(Rb, r_wm(end,:));

    kern4 = zeros([GRID_SNR 5]);
    for j = 1:5, kern4(:,:,:,j) = K_WM(j); end

    % ---------------- SMI settings ----------------
    SET = {};
    SET{end+1} = struct('name','SMI base', 'reg', SMI_BASE);
    for io = 1:numel(SMI_OAT)
        fld = SMI_OAT{io}{1};
        for v = SMI_OAT{io}{2}
            rg = SMI_BASE; rg.(fld) = v;
            SET{end+1} = struct('name', sprintf('SMI %s=%g', fld, v), 'reg', rg);
        end
    end
    for ie = 1:numel(SMI_EXTRA)
        rg = SMI_BASE; fn = fieldnames(SMI_EXTRA{ie}); nm = 'SMI';
        for k = 1:numel(fn)
            rg.(fn{k}) = SMI_EXTRA{ie}.(fn{k});
            nm = [nm sprintf(' %s=%g', fn{k}, SMI_EXTRA{ie}.(fn{k}))];
        end
        SET{end+1} = struct('name', nm, 'reg', rg);
    end
    gA = SMI_GRID{1}; vA = SMI_GRID{2}; gB = SMI_GRID{3}; vB = SMI_GRID{4};
    for a = vA
        for bq = vB
            rg = SMI_BASE; rg.(gA) = a; rg.(gB) = bq;
            SET{end+1} = struct('name', sprintf('SMI %s=%g %s=%g', gA,a,gB,bq), 'reg', rg);
        end
    end

    for k = 1:numel(SET)
        SH = zeros(NVOX, nc); t0 = tic;
        for is = 1:NSNR
            rows = find(snr_id == is);
            dwin = reshape(S_noisy(rows,:), [GRID_SNR Ndwi]);
            plm = SMI.get_plm_from_S_and_kernel(dwin, [0 Lf Lf Lf], kern4, ...
                    true(GRID_SNR), bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                    bvecs, CS_PHASE, D_FW, SET{k}.reg);
            p = reshape(plm, [NVOX_SNR nc-1]);
            s = [ones(NVOX_SNR,1) p].*repmat(sc', NVOX_SNR, 1);
            s(~isfinite(s)) = 0;
            SH(rows,:) = s;
        end
        JOB{end+1} = struct('arm','SMI','name',SET{k}.name,'Lmax',Lf,'SH',SH,'secs',toc(t0));
        fprintf('   Lmax %d  %-38s %6.1f s\n', Lf, SET{k}.name, toc(t0));
    end

    % ---------------- SSST-CSD settings ----------------
    for nl = CSD_NEG
      for ml = CSD_NORM
        for th = CSD_THR
          for ifl = 1:numel(CSD_FILT)
            ex = sprintf('-neg_lambda %g -norm_lambda %g -threshold %g', nl, ml, th);
            nm = sprintf('SSST neg=%g norm=%g thr=%g', nl, ml, th);
            if ~isempty(CSD_FILT{ifl})
                ex = [ex ' -filter ' CSD_FILT{ifl}];
                nm = [nm ' filt=' CSD_FILT{ifl}];
            end
            fo = fullfile(WDIR,'csdfod'); t0 = tic;
            [st,tx] = system(sprintf(['dwi2fod csd "%s.mif" "%s" "%s.mif" -mask "%s.mif" ' ...
                                      '%s -force -quiet 2>&1'], f_b3, Rb, fo, f_mask, ex));
            if st ~= 0
                fprintf('   SKIP %-38s %s\n', nm, strtrim(tx(1:min(60,end)))); continue
            end
            V = MR.read([fo '.mif']);
            % Two statements, not reshape(...)(:,1:nc): indexing directly into a
            % function result is Octave-only and is a syntax error in MATLAB.
            Vm = reshape(V, [NVOX size(V,4)]);
            JOB{end+1} = struct('arm','SSST-CSD','name',nm,'Lmax',Lf, ...
                                'SH',Vm(:,1:nc),'secs',toc(t0));
          end
        end
      end
    end

    % ---------------- MSMT-CSD settings ----------------
    for nl = MSMT_NEG
      for ml = MSMT_NORM
        nm = sprintf('MSMT neg=%g norm=%g', nl, ml); t0 = tic;
        fw_=fullfile(WDIR,'mw'); fg_=fullfile(WDIR,'mg'); fk_=fullfile(WDIR,'mc');
        [st,tx] = system(sprintf(['dwi2fod msmt_csd "%s.mif" "%s" "%s.mif" "%s" "%s.mif" ' ...
            '"%s" "%s.mif" -mask "%s.mif" -neg_lambda %g -norm_lambda %g -force -quiet 2>&1'], ...
            f_dwi, Rw, fw_, Rg, fg_, Rc, fk_, f_mask, nl, ml));
        if st ~= 0
            fprintf('   SKIP %-38s %s\n', nm, strtrim(tx(1:min(60,end)))); continue
        end
        V = MR.read([fw_ '.mif']);
        Vm = reshape(V, [NVOX size(V,4)]);          % see the note above
        JOB{end+1} = struct('arm','MSMT-CSD','name',nm,'Lmax',Lf, ...
                            'SH',Vm(:,1:nc),'secs',toc(t0));
      end
    end
    fprintf('   Lmax %d: MRtrix settings done\n', Lf);
end
fprintf('   %d settings collected\n\n', numel(JOB));

%% Step 3 -- score every collected fODF through one identical path
% Peaks come from sh2peaks, as in notebooks/smi_wm_60deg.m, so a number here is
% comparable with a number there. Shape metrics are against the band-limited
% truth at the SAME Lmax, so the band limit is not charged as method error.

ERR=nan(numel(JOB),NSNR); COR=nan(numel(JOB),NSNR); SPU=nan(numel(JOB),NSNR);
ACC=nan(numel(JOB),NSNR); REL=nan(numel(JOB),NSNR); NEG=nan(numel(JOB),NSNR);

fprintf('Step 3: scoring\n');
for k = 1:numel(JOB)
    Lf = JOB{k}.Lmax; nc = (Lf/2+1)*(Lf+1);
    SH = JOB{k}.SH;
    Ye = SMI.get_even_SH(de, Lf, CS_PHASE);
    TR = sh_gt(:,1:nc);

    % ---- peaks
    fo = fullfile(WDIR,'score');
    MR.write([fo '.mif'], reshape(SH,[GRID_ALL nc]));
    [st,tx] = system(sprintf('sh2peaks "%s.mif" "%s_p.mif" -num %d -force -quiet 2>&1', ...
                             fo, fo, PEAK_NUM));
    if st ~= 0, fprintf(2,'%s\n',tx); error('sh2peaks failed on %s', JOB{k}.name); end
    P  = reshape(MR.read([fo '_p.mif']), [NVOX 3*PEAK_NUM]);
    a0 = SH(:,1)/sqrt(4*pi);
    amp = zeros(NVOX,PEAK_NUM);
    for q = 1:PEAK_NUM, amp(:,q) = sqrt(sum(P(:,(q-1)*3+(1:3)).^2,2)); end
    ani = amp - repmat(a0,1,PEAK_NUM); ani(~isfinite(ani)) = -Inf;
    keep = (ani >= PEAK_REL*repmat(max(ani,[],2),1,PEAK_NUM)) & isfinite(amp) & (amp>0);
    nf = sum(keep,2); ae = nan(NVOX,1);
    for q = 1:NVOX
        if ~keep(q,1), continue, end
        v = P(q,1:3); nv = norm(v);
        if nv == 0 || ~isfinite(nv), continue, end
        v = v/nv; r = orient_id(q);
        ae(q) = min([acosd(min(abs(v*AX{r,1}(:)),1)), acosd(min(abs(v*AX{r,2}(:)),1))]);
    end

    % ---- shape, both sides put on p_00 = 1 so the comparison is scale free
    good = SH(:,1) > 0;
    SHn  = SH;
    SHn(good,:) = SH(good,:) .* repmat((1/sqrt(4*pi))./SH(good,1), 1, nc);
    SHn(~good,:) = NaN;
    A  = SHn*Ye';
    TA = TR*Ye';
    accv = nan(NVOX,1); relv = nan(NVOX,1); negv = nan(NVOX,1);
    for q = 1:NVOX
        if ~good(q), continue, end
        r = orient_id(q);
        u = SHn(q,2:end); t = TR(r,2:end);
        accv(q) = (u*t')/(norm(u)*norm(t));
        relv(q) = norm(A(q,:)-TA(r,:))/norm(TA(r,:));
        negv(q) = mean(A(q,:) < 0);
    end

    for is = 1:NSNR
        sel = (snr_id == is);
        ERR(k,is) = mean(ae(sel & isfinite(ae)));
        COR(k,is) = 100*mean(nf(sel) == 2);
        SPU(k,is) = mean(max(nf(sel)-2,0));
        ACC(k,is) = mean(accv(sel & isfinite(accv)));
        REL(k,is) = mean(relv(sel & isfinite(relv)));
        NEG(k,is) = mean(negv(sel & isfinite(negv)));
    end
end
fprintf('   done\n');

%% Step 4 -- the tables
fprintf('\n================ RESULTS ================\n');
fprintf('PEAK : err mean angular error (deg), cor%% true fibre count, spur per voxel\n');
fprintf('SHAPE: ACC angular correlation l>=2 (1 perfect), L2 relative error,\n');
fprintf('       neg fraction of the sphere below zero\n\n');
arms = {'SMI','SSST-CSD','MSMT-CSD'};
for iL = 1:numel(LMAX_LIST)
    for is = 1:NSNR
        fprintf('--- Lmax %d, SNR %s ---\n', LMAX_LIST(iL), SNRLAB{is});
        fprintf('  %-38s %7s %7s %7s | %7s %7s %7s\n', ...
                'setting','err','cor%','spur','ACC','L2','neg');
        for ia = 1:numel(arms)
            for k = 1:numel(JOB)
                if ~strcmp(JOB{k}.arm,arms{ia}) || JOB{k}.Lmax ~= LMAX_LIST(iL), continue, end
                fprintf('  %-38s %7.2f %7.1f %7.3f | %7.4f %7.4f %7.3f\n', ...
                        JOB{k}.name, ERR(k,is), COR(k,is), SPU(k,is), ...
                        ACC(k,is), REL(k,is), NEG(k,is));
            end
        end
        fprintf('\n');
    end
end

%% Step 5 -- winners, and where the two metric families disagree
fprintf('================ BEST SETTING PER METRIC ================\n');
fprintf('Where the left and right halves pick different settings, the peak and\n');
fprintf('shape families disagree and the choice depends on what consumes the fODF.\n\n');
for iL = 1:numel(LMAX_LIST)
  for is = 1:NSNR
    for ia = 1:numel(arms)
      idx = [];
      for k = 1:numel(JOB)
        if strcmp(JOB{k}.arm,arms{ia}) && JOB{k}.Lmax == LMAX_LIST(iL), idx(end+1) = k; end
      end
      if isempty(idx), continue, end
      [ve,ie] = min(ERR(idx,is)); [vc,ic] = max(COR(idx,is));
      [va,ja] = max(ACC(idx,is)); [vl,jl] = min(REL(idx,is));
      fprintf('Lmax %d SNR %-4s %-10s\n', LMAX_LIST(iL), SNRLAB{is}, arms{ia});
      fprintf('   best angular error : %-38s %.2f deg\n', JOB{idx(ie)}.name, ve);
      fprintf('   best correct count : %-38s %.1f%%\n',   JOB{idx(ic)}.name, vc);
      fprintf('   best ACC           : %-38s %.4f\n',     JOB{idx(ja)}.name, va);
      fprintf('   best relative L2   : %-38s %.4f\n',     JOB{idx(jl)}.name, vl);
      if ~strcmp(JOB{idx(ie)}.name, JOB{idx(jl)}.name)
          fprintf('   ** PEAK and SHAPE metrics disagree **\n');
      end
      fprintf('\n');
    end
  end
end

if ~KEEP_FILES && exist(WDIR,'dir')
    d = dir(fullfile(WDIR,'*'));
    for k = 1:numel(d)
        if ~d(k).isdir, delete(fullfile(WDIR,d(k).name)); end
    end
end
fprintf('=== sweep done ===\n');
