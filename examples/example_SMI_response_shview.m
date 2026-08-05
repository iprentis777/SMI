% examples/example_SMI_response_shview.m
%
% Viewing the SMI response kernel as zonal harmonics, the way MRtrix3's
% `shview` displays a response function.
%
% WHY THIS EXISTS
%
% SMI does not have a "response function". It has a kernel: a parametric
% Standard Model compartment description [f Da Depar Deperp fw] fitted per
% voxel, from which the rotational invariants K_l(b) follow analytically. CSD
% has the opposite arrangement -- a non-parametric response function, estimated
% once for the whole brain and stored as zonal harmonic coefficients in a text
% file.
%
% They describe the same object. For a single fibre pointing along z the SMI
% forward model (SMI.m:818-820) gives
%
%     R(theta) = sum_l K_l(b) (2l+1) P_l(cos theta) = sum_l r_l Y_l0(theta)
%     r_l      = K_l(b) sqrt((2l+1)*4*pi)
%
% and r_l is exactly what MRtrix keeps in a response function file. So an SMI
% kernel can be written out as an MRtrix response and handed to `shview`, and
% an MRtrix response can be read back and compared against an SMI kernel in the
% same units. This script does both, and draws the glyphs itself so the
% comparison can be made without leaving MATLAB.
%
% WHAT IT DRAWS
%
%   figure 1  the kernel four ways: K_l(b), the zonal coefficients per shell,
%             the angular profile R(theta), and the shview glyph per shell
%   figure 2  the compartments taken apart -- what intra-axonal, extra-axonal
%             and free water each contribute to the response
%   figure 3  a response glyph next to an fODF glyph, both from the same
%             renderer, so the deconvolution can be read off the picture
%
% and it writes `response_SMI_wm.txt`, which `shview response_SMI_wm.txt` opens.
%
% Helper functions live in helpers/SMI_response_helpers.m rather than at the end of
% this script, so that it runs under Octave as well as MATLAB.

clear; close all

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);
addpath(fullfile(repo_root, 'helpers'));

H  = SMI_response_helpers();
HF = fODF_modulation_helpers();

% ------------------------------------------------------------------ inputs
% A typical healthy white matter kernel [f Da Depar Deperp fw], diffusivities
% in um^2/ms. Replace with your own fit:
%
%     out    = SMI.fit(dwi, options);
%     kernel = H.kernel_from_out(out, [x y z]);
%
kernel   = [0.60 2.0 2.0 0.50 0.05];
bvals    = [0 1 2 3];        % ms/um^2, the shells to draw
Lmax     = 8;                % the SM kernel carries Legendre terms to l = 8
D_FW     = 3;                % um^2/ms
CS_phase = 1;                % SMI.fit's default
S0       = 1;                % scale of the written response file
RESPONSE_FILE = '';          % optional MRtrix response .txt to overlay

l_list = 0:2:Lmax;
Nsh    = numel(bvals);
lbl    = cell(1,Nsh);
Klbl   = cell(1,numel(l_list));
rlbl   = cell(1,numel(l_list));
llbl   = cell(1,numel(l_list));
for i = 1:Nsh, lbl{i} = sprintf('b = %g', bvals(i)); end
for j = 1:numel(l_list)
    Klbl{j} = sprintf('K_%d', l_list(j));
    rlbl{j} = sprintf('r_%d', l_list(j));
    llbl{j} = sprintf('l = %d', l_list(j));
end

% -------------------------------------------------- kernel -> zonal harmonics
[R_zh, K_l] = H.zh(kernel, bvals, Lmax, D_FW, S0);

fprintf('\nSMI kernel  [f Da Depar Deperp fw] = %s\n', mat2str(kernel,3));
fprintf('\nrotational invariants K_l(b)\n');
fprintf('%8s', 'b'); fprintf('%12s', Klbl{:}); fprintf('\n');
for i = 1:Nsh
    fprintf('%8g', bvals(i)); fprintf('%12.5f', K_l(i,:)); fprintf('\n');
end
fprintf('\nzonal harmonic coefficients r_l = K_l*sqrt((2l+1)*4pi)');
fprintf('   (the rows of an MRtrix response file)\n');
fprintf('%8s', 'b'); fprintf('%12s', rlbl{:}); fprintf('\n');
for i = 1:Nsh
    fprintf('%8g', bvals(i)); fprintf('%12.5f', R_zh(i,:)); fprintf('\n');
end

% ------------------------------------------------------------ self check
% The glyph is only worth looking at if it is the same object the fit uses.
% Synthesise the signal of a delta fODF along z with SMI's own forward model
% and compare it against the zonal reconstruction at the same directions.
%
% p_lm of a delta along z: only m = 0 survives and equals 1 in the normalized
% convention. In SMI's ordering the l block runs m = -l..l, so the m = 0 entry
% of block l sits at offset l inside that block.
plm_delta = zeros(1, Lmax*(Lmax+3)/2);
ofs = 0;
for il = 2:2:Lmax
    plm_delta(ofs + il + 1) = 1;
    ofs = ofs + 2*il + 1;
end

dirs_chk = HF.dirs(200);
max_err  = 0;
for i = 1:Nsh
    b_i  = bvals(i)*ones(1,size(dirs_chk,1));
    S_fw = HF.signal(plm_delta, [kernel 1 1], b_i, ones(size(b_i)), ...
                     zeros(size(b_i)), dirs_chk, Lmax, CS_phase, D_FW);
    S_zh = H.profile(R_zh(i,:)/S0, acos(min(1,max(-1,dirs_chk(:,3)))));
    max_err = max(max_err, max(abs(S_fw(:) - S_zh(:))));
end
fprintf('\nself check: zonal reconstruction vs SMI forward model, max|err| = %.2e\n', max_err);
if max_err > 1e-10
    error('the zonal harmonic response does not reproduce the SMI forward model')
end

% ---------------------------------------------------- MRtrix response file
H.write_response('response_SMI_wm.txt', R_zh);
fprintf('wrote response_SMI_wm.txt  --  view it with:  shview response_SMI_wm.txt\n');

R_ext = [];
if ~isempty(RESPONSE_FILE) && exist(RESPONSE_FILE,'file')
    R_ext = H.read_response(RESPONSE_FILE);
    fprintf('read %s : %d shells x %d coefficients\n', ...
            RESPONSE_FILE, size(R_ext,1), size(R_ext,2));
end

% ===================================================================
% FIGURE 1 -- the kernel, four ways
% ===================================================================
theta = linspace(0, pi, 361);
cols  = lines(max([Nsh numel(l_list) 4]));

figure('Name','SMI response kernel','Position',[80 80 1180 760]);

subplot(2,2,1); hold on
for j = 1:numel(l_list)
    plot(bvals, K_l(:,j), '-o', 'LineWidth', 1.6, 'Color', cols(j,:));
end
xlabel('b  (ms/\mum^2)'); ylabel('K_l(b)'); grid on
title('rotational invariants of the kernel')
legend(llbl, 'Location', 'northeast');

subplot(2,2,2);
bar(R_zh');
set(gca, 'XTick', 1:numel(l_list), 'XTickLabel', rlbl);
ylabel('r_l'); grid on
title('zonal harmonic coefficients (the MRtrix response rows)')
legend(lbl, 'Location', 'northeast');

subplot(2,2,3); hold on
for i = 1:Nsh
    plot(theta*180/pi, H.profile(R_zh(i,:), theta), 'LineWidth', 1.8, ...
         'Color', cols(i,:));
end
for i = 1:size(R_ext,1)
    plot(theta*180/pi, H.profile(R_ext(i,:), theta), '--', 'LineWidth', 1.2, ...
         'Color', cols(min(i,size(cols,1)),:));
end
xlabel('angle from the fibre  (deg)'); ylabel('R(\theta)'); grid on
xlim([0 180]); title('angular profile of the response')
legend(lbl, 'Location', 'northeast');

% The glyph row: one shview-style surface per shell, all to a common scale, so
% that the b dependence shows up as size the way shview shows it.
subplot(2,2,4); hold on
rmax = 0;
for i = 1:Nsh
    rmax = max(rmax, max(abs(H.profile(R_zh(i,:), theta))));
end
sep = 2.4;
for i = 1:Nsh
    [X,Y,Z,C] = H.zh_glyph(R_zh(i,:), 121, 121, 1/rmax);
    surf(X + (i-1)*sep, Y, Z, C, 'EdgeColor', 'none');
    text((i-1)*sep, 0, 1.35, lbl{i}, 'HorizontalAlignment', 'center');
end
axis equal off; view(0,0); camlight headlight; lighting gouraud
colormap(gca, 'jet'); title('response glyph per shell (as shview draws it)')

% ===================================================================
% FIGURE 2 -- what each compartment contributes
% ===================================================================
% The intra-axonal stick, the extra-axonal zeppelin and the free water pool are
% three response functions that add. Free water is isotropic, so it only ever
% moves r_0: it inflates the l = 0 coefficient and flattens the response
% without changing its shape at l >= 2. That is the reason SMI's fODF
% amplitude is blind to free water while an AFD-style fODF is not
% (REPORT_fODF_freewater.md).
f = kernel(1); Da = kernel(2); Depar = kernel(3); Deperp = kernel(4); fw = kernel(5);
comp_name   = {'intra-axonal', 'extra-axonal', 'free water'};
comp_kernel = {[1 Da 0 0 0], [0 Da Depar Deperp 0], [0 Da Da Da 1]};
comp_weight = [f, 1-f-fw, fw];

figure('Name','SMI kernel by compartment','Position',[100 100 1180 420]);
b_show = bvals(end);
for ic = 1:3
    subplot(1,4,ic);
    r_c = comp_weight(ic)*H.zh(comp_kernel{ic}, b_show, Lmax, D_FW, S0);
    plot(theta*180/pi, H.profile(r_c, theta), 'LineWidth', 1.8);
    xlabel('angle (deg)'); ylabel('R(\theta)'); grid on; xlim([0 180])
    title(sprintf('%s (weight %.2f)', comp_name{ic}, comp_weight(ic)))
end
subplot(1,4,4); hold on
r_tot = H.zh(kernel, b_show, Lmax, D_FW, S0);
plot(theta*180/pi, H.profile(r_tot, theta), 'k', 'LineWidth', 2.2);
for ic = 1:3
    r_c = comp_weight(ic)*H.zh(comp_kernel{ic}, b_show, Lmax, D_FW, S0);
    plot(theta*180/pi, H.profile(r_c, theta), '--', 'LineWidth', 1.4);
end
xlabel('angle (deg)'); ylabel('R(\theta)'); grid on; xlim([0 180])
title(sprintf('sum at b = %g', b_show))
legend({'total','intra','extra','free water'}, 'Location', 'northeast');

% ===================================================================
% FIGURE 3 -- response next to fODF, same renderer
% ===================================================================
% A 60 degree crossing, band limited to Lmax 6, drawn through the same code
% path as the response above. This is the pair `shview` shows when a response
% and an FOD are opened side by side.
Lmax_fod = 6;
dq  = HF.dirs(2000);
n1  = [0 0 1];
n2  = [sind(60) 0 cosd(60)];
fod = HF.watson(dq, n1, 16) + HF.watson(dq, n2, 16);
plm_cross = HF.mixture_plm(fod, dq, Lmax_fod, CS_phase);

figure('Name','response and fODF','Position',[120 120 900 430]);
subplot(1,2,1);
[X,Y,Z,C] = H.zh_glyph(R_zh(end,:), 121, 121, 1/rmax);
surf(X,Y,Z,C,'EdgeColor','none'); axis equal off; view(0,0)
camlight headlight; lighting gouraud; colormap(gca,'jet')
title(sprintf('response, b = %g', bvals(end)))

subplot(1,2,2);
[X,Y,Z,C] = H.sh_glyph(plm_cross, Lmax_fod, CS_phase, 91, 121, 1);
surf(X,Y,Z,C,'EdgeColor','none'); axis equal off; view(0,0)
camlight headlight; lighting gouraud; colormap(gca,'jet')
title('fODF, 60 deg crossing, Lmax 6')

fprintf('\ndone.\n');
