% dump_reference.m -- write the Octave side of the Python/Octave comparison.
%
% Every array this file writes is produced by the SAME functions the manuscript
% simulation uses -- SMI.get_even_SH, fODF_modulation_helpers, mc_config's
% protocol loader -- at the manuscript's own configuration. check_python_vs_octave.py
% recomputes each one in Python and reports max|err|.
%
% This is what makes "the Python notebook simulates the same experiment as the
% .m file" a measurement rather than a claim. It is deliberately NOT a test of
% the noise: the two languages draw from different generators, and the noise-free
% signal is the thing that has to be identical.
%
%   octave-cli --no-gui -q dump_reference.m
%
% Writes into data/ (gitignored), which is where binio.m and smi_sim.bin_load
% both look.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'oct_path.m'));

MC = mc_config();
H  = fODF_modulation_helpers();
RH = SMI_response_helpers();
IO = binio();

% ---- the manuscript configuration, copied from notebooks/smi_manuscript_60deg.m
K_WM     = [0.60 2.0 2.0 0.50 0.02];
D_FW     = 3;
KAPPA    = 16;
ANGLE    = 60;
AXIS1    = [0.30 -0.50 0.81]; AXIS1 = AXIS1/norm(AXIS1);
LMAX_GT  = 8;
CS_PHASE = 0;
PROTOCOL = 'hcp_real_3shell.txt';
B0_SNAP  = 0.05;
NDIR_Q   = 3000;
NDIR_E   = 1500;

fprintf('=== dumping the Octave reference for the Python comparison ===\n');

% ---- Step 1: the protocol, exactly as every arm reads it
[bvals, bvecs] = MC.load_protocol_file(PROTOCOL);
bvals(bvals < B0_SNAP) = 0;
IO.save('ref_bvals', bvals(:));
IO.save('ref_bvecs', bvecs);
fprintf('  protocol      %d volumes, %d distinct b\n', numel(bvals), numel(unique(bvals)));

% ---- the two direction sets
dq = H.dirs(NDIR_Q);
de = H.dirs(NDIR_E);
IO.save('ref_dq', dq);
IO.save('ref_de', de);

% ---- Step 2: the SH basis itself, on both direction sets and on the gradients
IO.save('ref_Y_dq_L8', SMI.get_even_SH(dq,    LMAX_GT, CS_PHASE));
IO.save('ref_Y_de_L6', SMI.get_even_SH(de,    6,       CS_PHASE));
IO.save('ref_Y_g_L8',  SMI.get_even_SH(bvecs, LMAX_GT, CS_PHASE));
% and at CS_phase = 1, because the two conventions are where sign errors hide
IO.save('ref_Y_de_L6_cs1', SMI.get_even_SH(de, 6, 1));
fprintf('  SH basis      Lmax %d on %d and %d directions, both CS conventions\n', ...
        LMAX_GT, NDIR_Q, NDIR_E);

% ---- Step 2: the ground truth fODF
axes_gt = {AXIS1, MC.rotate_about(AXIS1, ANGLE)};
fodf = zeros(size(dq,1), 1);
for k = 1:numel(axes_gt)
    fodf = fodf + H.watson(dq, axes_gt{k}, KAPPA);
end
plm_gt = H.mixture_plm(fodf, dq, LMAX_GT, CS_PHASE);
IO.save('ref_axis2',  axes_gt{2}(:));
IO.save('ref_fodf_q', fodf(:));
IO.save('ref_plm_gt', plm_gt(:));
fprintf('  ground truth  %d deg crossing, kappa %g, %d plm coefficients\n', ...
        ANGLE, KAPPA, numel(plm_gt));

% ---- the Watson zonal coefficients, which are what makes the CSD response
% dispersion matched rather than a delta
w1   = H.watson(dq, [0 0 1], KAPPA);
plmw = H.mixture_plm(w1, dq, LMAX_GT, CS_PHASE);
Lw   = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';
Mw   = [];
for il = 0:2:LMAX_GT, Mw = [Mw, -il:il]; end %#ok<AGROW>
shw  = [1; plmw(:)];
IO.save('ref_watson_pl', shw(Mw(:) == 0));
fprintf('  watson p_l    kappa %g, l = 0..%d\n', KAPPA, LMAX_GT);

% ---- Step 3: the kernel invariants and both responses
Kmat = RH.Kell(K_WM, bvals, LMAX_GT, D_FW);
IO.save('ref_Kell', Kmat);
IO.save('ref_resp_delta', RH.zh(K_WM, bvals, LMAX_GT, D_FW));
fprintf('  kernel        K_l over %d b values, l = 0..%d\n', numel(bvals), LMAX_GT);

% ---- Step 4: the noise-free signal, both ways
S_clean = H.signal(plm_gt, [K_WM 1 1], bvals, ones(1,numel(bvals)), ...
                   zeros(1,numel(bvals)), bvecs, LMAX_GT, CS_PHASE, D_FW);
IO.save('ref_S_clean', S_clean(:));

f_ = K_WM(1); Da_ = K_WM(2); Dep_ = K_WM(3); Dpp_ = K_WM(4); fw_ = K_WM(5);
cosang  = bvecs * dq';
Kdir    = f_         * exp(-bvals(:) .* (Da_*cosang.^2)) + ...
          (1-f_-fw_) * exp(-bvals(:) .* (Dpp_ + (Dep_-Dpp_)*cosang.^2)) + ...
          fw_        * exp(-bvals(:) * D_FW);
w = fodf/sum(fodf);
IO.save('ref_S_direct', (Kdir*w));
fprintf('  signal        noise free, harmonic and direct convolution\n');

% ---- pick_grid, whose factorisation the Python side must agree with
for N = [81 2187 3000 27]
    IO.save(sprintf('ref_grid_%d', N), MC.pick_grid(N)(:));
end
fprintf('  pick_grid     81, 2187, 3000, 27\n');

fprintf('=== done, %s ===\n', IO.dir());
