% =========================================================================
% Example: choosing the fODF deconvolution regularization parameters
%
% examples/example_fODF_regularization.m compares a few fixed settings. This script
% sweeps them, so that the values in the documentation are measured rather
% than guessed, and so that they can be re-measured for a different
% protocol by editing the protocol block below.
%
% The fODF is estimated from synthetic two-fiber voxels whose ground truth
% plm are known exactly, over a range of crossing angles and SNRs, and the
% estimate is scored against that ground truth. The sweep is staged:
%
%   1. lambda_nonneg
%   2. tau, at the best point of that sweep
%   3. Lmax_init, at the best point of that sweep
%
% Stages 1 and 4 were once a joint lambda_nonneg x lambda_tikhonov grid and a
% choice of Tikhonov matrix. Tikhonov damping has been removed from the
% toolbox as inert (Archive/patch_history/fODF_tikhonov/), so stage 1 is one
% dimensional and stage 4 is gone. Nothing else about the sweep changed --
% lambda_nonneg has NOT been re-optimized here.
%
% Scores (all averaged over crossing angles, SNRs and noise realizations):
%
%   rel err fODF   relative L2 error of the fODF over the sphere. This is
%                  the primary score: it is what "reconstructs the fODF
%                  similar to the ground truth" means, and it is
%                  dimensionless.
%   RMSE(plm)      error of the coefficients themselves
%   negative mass  fraction of the absolute mass of the fODF that is
%                  negative, i.e. how unphysical the estimate is
%   peak error     mean angle between each true fiber direction and the
%                  closest peak of the estimated fODF. This is the
%                  quantity tractography actually consumes.
%
% Scoring lives in helpers/fODF_regularization_score.m (a separate file rather
% than a local function, so that the script also runs in GNU Octave).
%
% All the deconvolutions for one parameter setting are batched into a
% single call to SMI.get_plm_from_S_and_kernel by stacking the noise
% realizations along the first image dimension, which is much faster than
% calling it once per realization and gives bit-identical results.
%
% Requires no data. Runtime is a few minutes.
%
% By: Santiago Coelho
% =========================================================================
clear; close all

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);
addpath(fullfile(repo_root, 'helpers'));

CS_phase = 1; D_FW = 3; Lmax = 6;
L_all = repelem(0:2:Lmax,2*(0:2:Lmax)+1);
L_rest = L_all(2:end);

% ---- Protocol: 3 shells with 64 directions each + 4 b0s ----------------
Ndirs_shell = 64;
dirs_shell = SMI.GetUniformHemisphereDirs(Ndirs_shell);
bvals = [0 1 2 3]; % [ms/um^2]
b = []; dirs = [];
for ii = 1:length(bvals)
    if bvals(ii)==0
        b = [b zeros(1,4)]; dirs = [dirs; repmat([0 0 1],4,1)];
    else
        b = [b bvals(ii)*ones(1,Ndirs_shell)]; dirs = [dirs; dirs_shell];
    end
end
beta = ones(size(b)); TE = zeros(size(b));
Ndwi = length(b); Lmax_shells = Lmax*ones(1,length(bvals));

% ---- Kernel ------------------------------------------------------------
f = 0.6; Da = 2; Depar = 2; Deperp = 0.5; fw = 0; T2a = 1; T2e = 1;
kernel_1vox = reshape([f Da Depar Deperp fw T2a T2e],[1 1 1 7]);

% ---- Conditions the parameters have to work across ---------------------
% A single crossing angle at a single SNR would give an optimum tuned to
% that one case. Averaging over a range gives settings that are not.
angles    = [40 60 90];   % crossing angle [degrees]
SNRs      = [20 30 50];
weights   = [0.6 0.4];    % volume fractions of the two fibers
smoothing = 0.05;         % makes the ground truth fODF non-negative
Nrep      = 30;           % noise realizations per condition

rng_seed = 1;
if exist('rng','builtin') || exist('rng','file'), rng(rng_seed), else, randn('seed',rng_seed); end

% ---- Build every voxel of the experiment once --------------------------
Ncond = numel(angles)*numel(SNRs);
Nvox  = Ncond*Nrep;
DWI   = zeros(Nvox,Ndwi);
PLM_T = zeros(Nvox,numel(L_rest));
COND  = zeros(Nvox,2);          % [angle index, SNR index]
FIB   = cell(numel(angles),1);
vv = 0;
for ia = 1:numel(angles)
    fibers = [1 0 0; cosd(angles(ia)) sind(angles(ia)) 0];
    FIB{ia} = fibers;
    plm_SH  = (weights*SMI.get_even_SH(fibers,Lmax,CS_phase)).*exp(-smoothing*L_all.*(L_all+1));
    plm_all = plm_SH./sqrt((2*L_all+1)/(4*pi));
    plm_a   = plm_all(2:end);
    S_a = SMI.SM_wFW_b_beta_TE_RealSphHarm(f,Da,Depar,Deperp,1-f-fw,T2a,T2e, ...
              plm_a,b,dirs,beta,TE,CS_phase,D_FW)';
    for is = 1:numel(SNRs)
        for rr = 1:Nrep
            vv = vv+1;
            DWI(vv,:)   = S_a + randn(1,Ndwi)/SNRs(is);
            PLM_T(vv,:) = plm_a;
            COND(vv,:)  = [ia is];
        end
    end
end
dwi4  = reshape(DWI,[Nvox 1 1 Ndwi]);
kernel = repmat(kernel_1vox,[Nvox 1 1 1]);
mask   = true(Nvox,1,1);
fprintf('Sweep over %d conditions x %d noise realizations = %d voxels per setting\n', ...
    Ncond, Nrep, Nvox);
fprintf('  crossing angles: %s deg | SNR: %s | Lmax = %d\n\n', ...
    mat2str(angles), mat2str(SNRs), Lmax);

% ---- Directions used to score the fODF ---------------------------------
dirs_odf = SMI.GetUniformHemisphereDirs(500);
Y_odf    = SMI.get_even_SH(dirs_odf,Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
AMP_T    = Y_odf(:,1) + Y_odf(:,2:end)*PLM_T.';        % [Ndirs x Nvox]
nrm_T    = sqrt(sum(AMP_T.^2,1));

% Neighbour list on the scoring directions, for peak detection
nNeighbors = 8;
G = abs(dirs_odf*dirs_odf.'); G(1:size(G,1)+1:end) = -Inf;
[~,ordG] = sort(G,2,'descend');
NB = ordG(:,1:nNeighbors);


% Peak amplitude of the ground truth, for reference on the figures
pk_true_mean = mean(max(AMP_T,[],1));
fprintf('  ground truth mean peak fODF amplitude: %.4f\n\n', pk_true_mean);

% =========================================================================
%  STAGE 1: lambda_nonneg
% =========================================================================
% This was a joint lambda_nonneg x lambda_tikhonov grid. Tikhonov damping has
% been removed from the toolbox as inert, so the sweep is one dimensional and
% lambda_nonneg is the only weight left to choose.
lam_nn  = [0 0.3 1 3 10 30 100]; % 0 = non-negativity constraint disabled

sz    = [1 numel(lam_nn)];
E_odf = nan(sz); E_plm = nan(sz); E_neg = nan(sz);
E_pk  = nan(sz); E_amp = nan(sz); E_rat = nan(sz); E_npk = nan(sz);

fprintf('STAGE 1: lambda_nonneg  (primary score: rel err fODF)\n');
fprintf('%12s%12s%12s\n','lambda_nonneg','rel err fODF','peak ratio');
for in = 1:numel(lam_nn)
    reg = struct();
    reg.flag_nonneg = double(lam_nn(in) > 0);
    if reg.flag_nonneg, reg.lambda_nonneg = lam_nn(in); end
    [E_odf(in),E_plm(in),E_neg(in),E_amp(in),E_rat(in),E_pk(in),E_npk(in)] = ...
        fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                      PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    fprintf('%12.3g%12.4f%12.3f\n',lam_nn(in),E_odf(in),E_rat(in));
end
fprintf('\n  peak amplitude ratio is peak(estimate)/peak(truth). 1 = no shrinkage\n');

[bestE,in_b] = min(E_odf);
fprintf('\n  best: lambda_nonneg = %g  -> rel err fODF %.4f\n', lam_nn(in_b), bestE);
fprintf('  at that point: RMSE(plm) %.4f | negative mass %.4f | peak error %.2f deg\n', ...
    E_plm(in_b), E_neg(in_b), E_pk(in_b));
fprintf('                 peak amplitude %.4f (%.1f%% of ground truth) | %.2f peaks/voxel\n\n', ...
    E_amp(in_b), 100*E_rat(in_b), E_npk(in_b));

best.flag_nonneg     = double(lam_nn(in_b) > 0);
best.lambda_nonneg   = max(lam_nn(in_b),1);

% =========================================================================
%  STAGE 2: tau
% =========================================================================
taus = [0.02 0.05 0.1 0.2 0.4];
E_tau = nan(size(taus)); R_tau = nan(size(taus));
fprintf('STAGE 2: tau (at the stage 1 optimum)\n');
for k = 1:numel(taus)
    reg = best; reg.tau = taus(k);
    [E_tau(k),~,~,~,R_tau(k)] = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                             PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    fprintf('  tau = %5.2f -> rel err fODF %.4f | peak ratio %.3f\n', taus(k), E_tau(k), R_tau(k));
end
[~,k_b] = min(E_tau); best.tau = taus(k_b);
fprintf('  best tau = %g\n\n', best.tau);

% =========================================================================
%  STAGE 3: Lmax_init
% =========================================================================
Lin = 2:2:Lmax;
E_lin = nan(size(Lin)); R_lin = nan(size(Lin));
fprintf('STAGE 3: Lmax_init (at the stage 1-2 optimum)\n');
for k = 1:numel(Lin)
    reg = best; reg.Lmax_init = Lin(k);
    [E_lin(k),~,~,~,R_lin(k)] = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                             PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    fprintf('  Lmax_init = %d -> rel err fODF %.4f | peak ratio %.3f\n', Lin(k), E_lin(k), R_lin(k));
end
[~,k_b] = min(E_lin); best.Lmax_init = Lin(k_b);
fprintf('  best Lmax_init = %d\n\n', best.Lmax_init);

% =========================================================================
%  Per SNR
% =========================================================================
fprintf('\nPer-SNR optimum over the stage 1 grid\n');
E_snr = nan(numel(SNRs),numel(lam_nn));
R_snr = nan(size(E_snr));
for is = 1:numel(SNRs)
    sel = COND(:,2)==is;
    for in = 1:numel(lam_nn)
        reg = struct();
        reg.flag_nonneg = double(lam_nn(in) > 0);
        if reg.flag_nonneg, reg.lambda_nonneg = lam_nn(in); end
        [E_snr(is,in),~,~,~,R_snr(is,in)] = ...
            fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                                  PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB,sel);
    end
    [e,ii] = min(E_snr(is,:));
    fprintf('  SNR %2d: lambda_nonneg = %-4g -> %.4f (peak ratio %.3f)\n', ...
        SNRs(is), lam_nn(ii), e, R_snr(is,ii));
end

% =========================================================================
%  Recommendation
% =========================================================================
fprintf('\n==================== RECOMMENDED ====================\n');
fprintf('options.fODF_regularization.flag_nonneg      = %d;\n', best.flag_nonneg);
fprintf('options.fODF_regularization.lambda_nonneg    = %g;\n', best.lambda_nonneg);
fprintf('options.fODF_regularization.tau              = %g;\n', best.tau);
fprintf('options.fODF_regularization.Lmax_init        = %d;\n', best.Lmax_init);
fprintf('=====================================================\n');

% =========================================================================
%  FIGURES
%
%  Written for a manuscript: panel letters, consistent styling, and each
%  figure saved to disk in a vector format alongside a raster preview.
%  Set saveFigures = false to only display them.
% =========================================================================
saveFigures = true;
figDir      = fullfile(pwd,'figures_fODF_sweep');
figFormats  = {'-dpng','-depsc'};   % raster preview + vector for the manuscript
figRes      = '-r300';
if saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10)
yt = 1:numel(lam_nn);
lamNNlab = arrayfun(@(v) num2str(v),lam_nn,'UniformOutput',false); lamNNlab{1} = 'off';
panel = @(s) text(-0.16,1.06,s,'Units','normalized','FontWeight','bold','FontSize',12);

% ---- Figure 1: the regularization curves, all four scores --------------
f1 = figure('Color','w','Name','F1 regularization curves','Position',[60 60 1000 760]);
maps = {E_odf,'relative fODF error','a'; E_neg,'negative mass','b'; ...
        E_pk,'peak angular error [deg]','c'; E_rat,'peak amplitude ratio','d'};
for k = 1:4
    subplot(2,2,k), hold on
    plot(1:numel(lam_nn),maps{k,1},'o-','LineWidth',1.4)
    plot(in_b,maps{k,1}(in_b),'r*','MarkerSize',13,'LineWidth',1.5)
    set(gca,'XTick',yt,'XTickLabel',lamNNlab), box on, grid on
    xlabel('\lambda_{nonneg}')
    title(maps{k,2}), panel(maps{k,3})
end

% ---- Figure 2: does more regularization shrink the peaks? --------------
f2 = figure('Color','w','Name','F2 peak amplitude shrinkage','Position',[60 60 900 400]);
subplot(1,2,1), hold on
plot(1:numel(lam_nn),E_rat,'o-','LineWidth',1.3)
plot([1 numel(lam_nn)],[1 1],'k--','HandleVisibility','off')
set(gca,'XTick',yt,'XTickLabel',lamNNlab), box on, grid on
xlabel('\lambda_{nonneg}'), ylabel('peak(estimate) / peak(truth)')
title('shrinkage vs non-negativity weight'), panel('a')

subplot(1,2,2), hold on
scatter(E_rat(:),E_odf(:),28,'filled')
plot(E_rat(in_b),E_odf(in_b),'r*','MarkerSize',14,'LineWidth',1.5)
yl = get(gca,'YLim'); plot([1 1],yl,'k--'), ylim(yl)
box on, grid on
xlabel('peak amplitude ratio'), ylabel('relative fODF error')
title('accuracy vs shrinkage trade-off'), panel('b')

% ---- Figure 3: the secondary parameters --------------------------------
f3 = figure('Color','w','Name','F3 secondary parameters','Position',[60 60 900 400]);
subplot(1,2,1), hold on
plot(taus,E_tau,'o-','LineWidth',1.4)
plot(best.tau,min(E_tau),'r*','MarkerSize',13,'LineWidth',1.5)
xlabel('\tau'), ylabel('relative fODF error'), box on, grid on
title('threshold \tau'), panel('a')

subplot(1,2,2), hold on
bar(Lin,E_lin,0.5)
xlabel('L_{max,init}'), ylabel('relative fODF error'), box on
set(gca,'XTick',Lin), title('initialization order'), panel('b')

% ---- Figure 4: SNR dependence ------------------------------------------
f4 = figure('Color','w','Name','F4 SNR dependence','Position',[60 60 560 400]);
hold on
for is = 1:numel(SNRs)
    plot(1:numel(lam_nn),E_snr(is,:),'o-','LineWidth',1.3, ...
         'DisplayName',sprintf('SNR %d',SNRs(is)))
end
set(gca,'XTick',yt,'XTickLabel',lamNNlab), box on, grid on
xlabel('\lambda_{nonneg}'), ylabel('relative fODF error')
legend('Location','best'), title('error vs non-negativity weight')

if saveFigures
    figs  = [f1 f2 f3 f4];
    names = {'F1_landscape','F2_peak_shrinkage','F3_secondary_params','F4_snr'};
    for k = 1:numel(figs)
        for ff = 1:numel(figFormats)
            print(figs(k), fullfile(figDir,names{k}), figFormats{ff}, figRes);
        end
    end
    fprintf('\nFigures written to %s\n', figDir);
end
