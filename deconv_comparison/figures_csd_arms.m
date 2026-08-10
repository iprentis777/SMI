function figures_csd_arms(outdir)
% figures_csd_arms -- manuscript figures over the CSD arms the Python notebook produced.
%
%   octave-cli --no-gui -q figures_csd_arms.m
%   >> figures_csd_arms                     % draw only
%   >> figures_csd_arms('../Figures')       % draw and print to PNG
%
% Reads everything notebooks/csd_manuscript_60deg.ipynb exported and draws it in
% the style of smi_manuscript_60deg.m's figures, so the CSD arms can go straight
% into the manuscript without re-running anything.
%
%   Figure 1  the response each arm was given, one glyph per shell
%   Figure 2  reconstructed fODFs: one row per arm, truth then one column per SNR
%   Figure 3  bias, spread and spurious peaks against SNR, one row per arm
%
% Every 3D panel opens in an ISOMETRIC view so all subfigures are readable
% without touching the camera, and glyphs use the shview convention the rest of
% this package uses: radius is |amplitude|, colour is the SIGNED amplitude, so
% negative lobes show as a colour change rather than being folded silently into
% the surface.
%
% ONLY FIGURE 1 USES A SHARED RADIAL SCALE, and getting that to work takes two
% things where only one is obvious: `axis equal` fixes a panel's aspect RATIO,
% not its LIMITS, so a glyph drawn at half the radius gets an axis range half as
% wide and lands on the page at exactly the same size. Figure 1 scales every
% glyph AND pins each panel. Figures 2 and 3 autoscale deliberately -- they are
% about shape, and shrinking one into illegibility would buy nothing.
%
% GRAPHICS ARE NOT GUARANTEED HERE. Under Octave in a headless container only
% the gnuplot toolkit is usually available, and it ignores camlight and
% lighting, so glyphs come out flat shaded; `print` may fail outright. Every
% figure is therefore wrapped, and a failure prints rather than taking the run
% down. The Python notebook's own figures ARE verified to render, and the
% printed tables below are drawn from the same arrays either way.

if nargin < 1, outdir = ''; end
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'oct_path.m'));

RH = SMI_response_helpers();
H  = fODF_modulation_helpers();

C = read_csd_export();
NCOND = numel(C.angles);
NSNR  = numel(C.snr_list);
NARM  = numel(C.arms);
NL    = numel(C.lmax_list);
[~, SNR_ORD] = sort(C.snr_list);           % ascending, Inf last

SNR_LABEL = cell(1, NSNR);
for i = 1:NSNR
    if isinf(C.snr_list(i)), SNR_LABEL{i} = 'inf';
    else, SNR_LABEL{i} = sprintf('%g', C.snr_list(i)); end
end

% ---- the printed tables, which are the deliverable whether or not graphics work
fprintf('\n=== CSD arms, as exported ===\n');
fprintf('kernel %s, kappa %g, %s response, protocol %s\n', ...
        mat2str(C.kernel), C.kappa, C.response_mode, C.protocol);
for ia = 1:NARM
    fprintf('\n  %s\n', C.arms{ia});
    for iL = 1:NL
        fprintf('    Lmax %d\n', C.lmax_list(iL));
        fprintf('      condition    SNR   correct count   mean err   std err   spurious\n');
        for ic = 1:NCOND
            for k = 1:NSNR
                is = SNR_ORD(k);
                fprintf('      %7d deg  %5s   %11.1f%%   %8.2f  %8.2f   %8.3f\n', ...
                        C.angles(ic), SNR_LABEL{is}, ...
                        C.scores(ia,iL,is,ic,1), C.scores(ia,iL,is,ic,2), ...
                        C.scores(ia,iL,is,ic,3), C.scores(ia,iL,is,ic,5));
            end
            fprintf('      %7d deg ceiling %10d peaks   %8.2f     (band limit at Lmax %d)\n', ...
                    C.angles(ic), C.ceiling(ia,iL,ic,1), C.ceiling(ia,iL,ic,2), ...
                    C.lmax_list(iL));
        end
    end
end
fprintf('\n');

% ---------------------------------------------------------------- Figure 1
try
    K   = C.kernel;
    b   = unique(round(C.bvals*1e6)/1e6);
    nsh = numel(b);
    % The response as the notebook built it: the kernel, times the Watson zonal
    % coefficients when the response is dispersion matched. Rebuilt here from
    % the manifest rather than exported, so this figure cannot silently show a
    % different response from the one the binaries were given.
    pl = watson_pl_local(H, C.kappa, C.lmax_gt, C.cs_phase);
    Rr = RH.zh(K, b, C.lmax_gt, 3);
    if strcmp(C.response_mode, 'dispersed')
        Rr = Rr .* repmat(pl(:)', size(Rr,1), 1);
    end

    % One shared radial scale: find the largest glyph radius first, scale every
    % glyph by it, THEN pin every panel to the same limits. Both steps are
    % needed; see the header.
    rmax = 0;
    for i = 1:nsh
        [X,Y,Z] = RH.zh_glyph(Rr(i,:), 61, 61);
        rmax = max(rmax, max(sqrt(X(:).^2 + Y(:).^2 + Z(:).^2)));
    end
    GL = 1.05;
    figure('Name','Figure 1 -- the response each CSD arm was given');
    for i = 1:nsh
        subplot(1, nsh, i);
        [X,Y,Z,Cc] = RH.zh_glyph(Rr(i,:), 61, 61, 1/rmax);
        surf(X,Y,Z,Cc,'EdgeColor','none'); axis equal off;
        set(gca,'XLim',[-GL GL],'YLim',[-GL GL],'ZLim',[-GL GL]);
        view(135,25);                                   % isometric
        title(sprintf('b = %.2f', b(i)));
    end
    annotate_fig(sprintf('Figure 1   %s response, one shared radial scale', ...
                         C.response_mode));
    printfig(outdir, 'csd_fig1_response');
catch err
    fprintf(2, '** FIGURE 1 FAILED ** %s\n', err.message);
end

% ---------------------------------------------------------------- Figure 2
try
    iL = 1; Lf = C.lmax_list(iL);
    nrow = NARM; ncol = 1 + NSNR;
    figure('Name','Figure 2 -- reconstructed fODFs, one row per arm');
    for ia = 1:NARM
        % Column 1: the ground truth, at the truth's own order.
        subplot(nrow, ncol, (ia-1)*ncol + 1);
        plm_t = sh_to_plm(C.truth_sh(1,:), C.lmax_gt);
        [X,Y,Z,Cc] = RH.sh_glyph(plm_t, C.lmax_gt, C.cs_phase, 61, 61);
        surf(X,Y,Z,Cc,'EdgeColor','none'); axis equal off; view(135,25);
        title(sprintf('%s: truth', C.arms{ia}), 'Interpreter','none');

        for k = 1:NSNR
            is   = SNR_ORD(k);
            rows = find(C.cond_id == 1 & C.snr_id == is);
            m    = mean(C.sh{ia}{iL}(rows,:), 1);
            % An MRtrix FOD is unnormalised, so it is put on the p_00 = 1
            % convention before the shared glyph renderer sees it. That is what
            % sh_to_plm does, and it preserves peak orientation exactly because
            % it scales every band by the same factor.
            subplot(nrow, ncol, (ia-1)*ncol + 1 + k);
            [X,Y,Z,Cc] = RH.sh_glyph(sh_to_plm(m, Lf), Lf, C.cs_phase, 61, 61);
            surf(X,Y,Z,Cc,'EdgeColor','none'); axis equal off; view(135,25);
            title(sprintf('SNR %s', SNR_LABEL{is}));
        end
    end
    annotate_fig(sprintf('Figure 2   %d deg crossing, Lmax %d, mean fODF per block', ...
                         C.angles(1), Lf));
    printfig(outdir, 'csd_fig2_fodf');
catch err
    fprintf(2, '** FIGURE 2 FAILED ** %s\n', err.message);
end

% ---------------------------------------------------------------- Figure 3
try
    fin = SNR_ORD(isfinite(C.snr_list(SNR_ORD)));
    xs  = C.snr_list(fin);
    mets = {2, 'mean angular error (deg)'; 3, 'std of angular error (deg)'; ...
            5, 'spurious peaks per voxel'};
    figure('Name','Figure 3 -- bias, spread and spurious peaks against SNR');
    for ia = 1:NARM
        for j = 1:3
            subplot(NARM, 3, (ia-1)*3 + j); hold on
            for iL = 1:NL
                ys = squeeze(C.scores(ia,iL,fin,1,mets{j,1}));
                plot(xs, ys(:), 'o-', 'LineWidth', 1.4, ...
                     'DisplayName', sprintf('Lmax %d', C.lmax_list(iL)));
            end
            set(gca,'XScale','log'); grid on
            xlabel('SNR'); ylabel(mets{j,2});
            if j == 1, title(C.arms{ia}, 'Interpreter','none'); end
            if j == 1 && ia == 1 && NL > 1, legend('show','Location','best'); end
        end
    end
    annotate_fig(sprintf('Figure 3   %d deg crossing, %s kernel, %s response', ...
                         C.angles(1), C.preset, C.response_mode));
    printfig(outdir, 'csd_fig3_summary');
catch err
    fprintf(2, '** FIGURE 3 FAILED ** %s\n', err.message);
end

fprintf('\nDone. The arrays behind every panel are C.scores and C.sh from\n');
fprintf('read_csd_export; test_csd_roundtrip.m checks them against the notebook.\n');
end

% =====================================================================
function L = deg_vec(Lmax)
% The l of every even-order SH coefficient, in SMI's ordering.
L = repelem(0:2:Lmax, 2*(0:2:Lmax)+1)';
end

% =====================================================================
function plm = sh_to_plm(sh, Lmax)
% plm = sh_to_plm(sh, Lmax)
%
% Put ANY set of even-order SH coefficients onto SMI's normalised convention:
% l = 2..Lmax only, with p_00 = 1 implied.
%
%   f_lm = plm .* sqrt((2l+1)/(4*pi)),   f_00 = 1/sqrt(4*pi)
%
% so rescaling sh by (1/sqrt(4*pi))/sh(1) makes its l = 0 term the convention's,
% and dividing out sqrt((2l+1)/(4*pi)) recovers plm. For an fODF that is already
% in the convention this is the identity; for an MRtrix FOD, whose l = 0 term is
% unnormalised and carries apparent fibre density, it is the change of scale
% that makes the two drawable side by side.
%
% Uniform scaling of every band preserves peak ORIENTATION exactly. (Clipping
% coefficients independently would not, and must never be done -- "README for
% Claude", section 3, item 4.)
sh   = sh(:);
L    = deg_vec(Lmax);
n    = min(numel(sh), numel(L));
sh   = sh(1:n); L = L(1:n);
scal = (1/sqrt(4*pi)) / sh(1);
plm  = (sh(2:end) * scal) ./ sqrt((2*L(2:end)+1)/(4*pi));
plm  = plm(:)';
end

% =====================================================================
function pl = watson_pl_local(H, kappa, Lmax, CS_phase)
% Normalised zonal coefficients p_l (p_0 = 1) of a Watson about +z, by the same
% sampled-grid projection the ground truth uses -- so the response drawn here
% and the response the binaries were handed cannot disagree about what a Watson
% at this kappa is.
dq  = H.dirs(3000);
plm = H.mixture_plm(H.watson(dq, [0 0 1], kappa), dq, Lmax, CS_phase);
sh  = [1; plm(:)];
m   = [];
for il = 0:2:Lmax, m = [m, -il:il]; end %#ok<AGROW>
pl  = sh(m(:) == 0);
end

% =====================================================================
function annotate_fig(txt)
% A figure-level title that works in both MATLAB and Octave. sgtitle exists only
% in recent MATLAB, so this falls back to an axes-free annotation.
try
    sgtitle(txt, 'Interpreter', 'none');
catch
    try
        % Kept to a thin strip at the very top: at [0 0.94 1 0.06] it collided
        % with the first row's own subplot titles under gnuplot.
        annotation('textbox', [0 0.965 1 0.035], 'String', txt, ...
                   'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
                   'Interpreter', 'none');
    catch
        fprintf('   %s\n', txt);
    end
end
end

% =====================================================================
function printfig(outdir, name)
% Print only when asked, and never take the run down if the toolkit cannot.
if isempty(outdir), return, end
if ~exist(outdir, 'dir'), mkdir(outdir); end
f = fullfile(outdir, [name '.png']);
try
    print(f, '-dpng', '-r150');
    fprintf('   wrote %s\n', f);
catch err
    fprintf(2, '   could not print %s: %s\n', f, err.message);
end
end
