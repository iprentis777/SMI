%% Manuscript figures: SMI on 60 degree crossings, healthy and edema
% The manuscript version of |smi_simulation_walkthrough.m|, cut down to the two
% configurations the manuscript uses -- *a single fibre and a 60 degree
% crossing* -- and swept across *a list of SNRs* rather than run at one.
% Everything else -- forward model, checks, fit, export -- is the same code.
% *Its configuration is local:* every knob lives in the Configuration block
% below, including the kernel, so retuning the experiment or switching to the
% edema kernel never means opening another file.
%
% *Why a sweep.* One SNR answers "does it work here". A sweep answers "where
% does it stop working", which is what a manuscript figure is for. It is also
% where the published comparison changes places: constrained SMI is the
% noise-robust arm and SSST-CSD the high-SNR angular-resolution one
% (|Reports/REPORT_SMI_deconvolution_MonteCarlo.md|, section 6.4). Nothing here
% assumes that ordering; the sweep is what would show it.
%
% *Why only 60 degrees.* 30 and 45 are kept in the general walkthrough. At
% |kappa = 16| a 30 degree crossing does not separate even in the noise-free
% ground truth, and 45 separates only from Lmax 6 upward, so neither is a clean
% figure. 60 separates at every Lmax and every kappa tested, which makes it the
% configuration where a difference between methods is attributable to the
% method. The single fibre condition is not simulated here either -- it is in
% the general walkthrough.
%
% *This file is meant to be read while it runs.* It is a plain |.m| script, so
% it runs as-is in MATLAB and in GNU Octave; MATLAB can also open it and save
% it as a Live Script (|.mlx|) with no edits, since the sections and markup
% below are publish-style.
%
% *Nothing here is a reimplementation.* The forward model comes from
% |helpers/fODF_modulation_helpers.m| and the fit is the real |SMI.fit|. The
% experiment's settings are in the Configuration block below rather than in
% |mc_config.m|, so this file can be retuned on its own; the geometry and
% protocol *utilities* are still shared with |gen_montecarlo.m| so those
% conventions cannot drift.
%
% *What you should end up believing.* Not "the code ran". Each step ends with
% one or more *CHECK* lines that compare its output against something computed
% a different way, and print |ok| or |** FAILED **| next to the number. If you
% read only the CHECK lines, you have audited the simulation.
%
%  Step 1  the acquisition protocol         -> a real HCP 3-shell scheme
%  Step 2  the ground truth fibre geometry  -> a single fibre and a 60 deg crossing
%  Step 3  the kernel, and it as a response -> K_l(b), zonal harmonics
%  Step 4  forward convolution              -> noise-free signal
%  Step 5  Rician noise, one block per SNR  -> the measured data
%  Step 6  SMI.fit at each Lmax and each SNR-> kernel and fODF recovery
%  Step 6b dwi2fod csd / msmt_csd           -> the two CSD arms, same signal
%  Step 7  peaks, angular error, spurious   -> where does it stop working?
%  Step 8  export every arm for MRtrix      -> check them outside this file
%  Fig 1   ground truth fibre geometry      -> the two configurations
%  Fig 2   the kernel as a response         -> one glyph per shell
%  Fig 3   the signal                       -> b down the rows, SNR across
%  Fig 4   SMI fODFs, healthy               -> truth, then SNR across
%  Fig 5   SMI fODFs, edema                 -> truth, then SNR across
%  Fig 6   bias, std and spurious vs SNR    -> one row per kernel AND arm
%  Fig 7   the arms side by side            -> one figure per kernel
%  Step 9  how to scale this up
%
% *Runtime is hours, not minutes, at the defaults below.* Step 6 is
% |numel(SNR_LIST) * numel(LMAX_LIST)| separate |SMI.fit| calls -- 18 at the
% defaults -- each on |2 * NREP| voxels. Two costs are now comparable rather
% than one dominating: |SMI.fit| trains its polynomial regression once per call
% regardless of voxel count, and at |NREP = 1000| the per-voxel constrained
% deconvolution matters too. *|NREP| is the knob*: drop it to 25 for a version
% that runs while you read it, and raise it once the output looks right. The
% file prints the elapsed time of every fit as it goes, so the first two lines
% of Step 6 tell you what the whole run will cost.

%% Configuration -- every knob in this simulation is here
% *Nothing below this block needs editing to retune the experiment*, and nothing
% in it has to be looked up elsewhere. The three things this file still takes
% from |mc_config.m| are stateless utilities, not settings: |pick_grid|,
% |rotate_about| and |load_protocol_file|. They are shared rather than copied so
% the fibre-axis convention and the protocol reader cannot drift away from
% |gen_montecarlo.m|.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
pkgdir = fileparts(here);                    % deconv_comparison/
if exist('OCTAVE_VERSION', 'builtin')        % quieten Octave's package-load noise
    warning('off', 'all');                   % so the CHECK lines are easy to find
end
run(fullfile(pkgdir, 'oct_path.m'));

MC   = mc_config();                          % utilities only, see above
H    = fODF_modulation_helpers();
RH   = SMI_response_helpers();
VERDICT = {'** FAILED **', 'ok'};            % VERDICT{1+condition}

% ---------------------------------------------------------------- the tissue
% BOTH kernels are simulated, back to back, in one run. Everything downstream
% carries the kernel name, so a healthy result and an edema result can never be
% mistaken for one another.
%
% SMI's compartment vector is [f Da Depar Deperp fw]:
%
%   f        intra-axonal (stick) fraction
%   Da       intra-axonal diffusivity along the stick        um^2/ms
%   Depar    extra-axonal diffusivity parallel to the fibre  um^2/ms
%   Deperp   extra-axonal diffusivity perpendicular          um^2/ms
%   fw       free water fraction, D fixed at D_FW below
%
% The extra-axonal fraction is whatever is left: 1 - f - fw.
%
% The edema kernel is shaped from a FITTED kernel rather than invented, which is
% why the diffusivities are not round numbers and why Depar exceeds Da -- a fit
% is free to land there and a synthetic kernel usually would not be written that
% way. Two deliberate departures from the numbers as fitted, both of which
% change what the result means:
%
% * *f is 0.10, not 0.05.* SMI's default training prior has its lower bound for
%   f at exactly 0.05 (SMI.m, MLTraining.bounds). A truth sitting ON the
%   boundary of the prior is estimated with a one-sided bias toward the prior
%   mean, which is indistinguishable from a real effect. 0.10 puts the truth
%   inside the prior so the result is about the data.
% * *fw is 0.35*, which the fitted kernel did not specify. Edema modelled as
%   added free water. The default prior caps fw at 0.5, so 0.35 is inside it.
%
% Deperp = 1.15 is still close to its own prior cap of 1.2. That one is left
% alone -- it is the number from the fit -- but if the extra-axonal
% diffusivities come back biased low, this is the first place to look.
KERNELS = { ...
    struct('name', 'healthy', 'K', [0.60 2.0 2.0 0.50 0.02], ...
           'note', 'the kernel the published Monte Carlo used'), ...
    struct('name', 'edema',   'K', [0.10 2.4 2.7 1.15 0.35], ...
           'note', 'low axonal fraction, high extra-axonal diffusivity, added free water') };
NKERN = numel(KERNELS);

D_FW  = 3;          % free water diffusivity, um^2/ms
KAPPA = 16;         % Watson concentration of each fibre population. Finite, not
                    % a delta: a response estimated from real white matter has
                    % already absorbed fibre dispersion, and a delta truth would
                    % create a mismatch that does not exist in practice.

% ------------------------------------------------------------- the geometry
ANGLES = 60;        % crossing angles in degrees. The manuscript uses the 60
                    % degree crossing only. Measured on the noise-free truth it
                    % separates at every Lmax and every KAPPA tested, where 30
                    % never separates at KAPPA = 16 and 45 only from Lmax 6 up.
                    % That is what makes a difference between methods
                    % attributable to the method rather than to the band limit.
                    % The single fibre condition is deliberately NOT simulated
                    % here; smi_simulation_walkthrough.m still covers it.
AXIS1  = [0.30 -0.50 0.81];       % first fibre axis, off every coordinate plane
AXIS1  = AXIS1/norm(AXIS1);

% ------------------------------------------------------------ the experiment
% *SMOKE_TEST is the knob to reach for first.* At the manuscript settings this
% file is NKERN x numel(SNR_LIST) x numel(LMAX_LIST) = 42 SMI.fit calls and runs
% in hours, which makes it useless for checking that a change still works. The
% smoke configuration is one Lmax and three SNRs and runs in minutes; every
% CHECK still executes, and the two MRtrix arms cost seconds either way because
% they scale with voxel count rather than with fit count.
%
% Its numbers are INDICATIVE ONLY -- 27 realisations is not a Monte Carlo -- and
% every printout says so.
SMOKE_TEST = false;

if SMOKE_TEST
    NREP      = 27;               % realisations per condition PER SNR
    SNR_LIST  = [10 30 Inf];
    LMAX_LIST = 6;
else
    NREP      = 1000;                     % realisations per condition PER SNR
    SNR_LIST  = [5 10 20 30 50 100 Inf];  % 1/sigma at each level, in any order
    LMAX_LIST = [4 6 8];                  % angular orders to fit at
end
LMAX_GT   = 8;                    % ground truth angular order. A CEILING, not a
                                  % choice: SMI's kernel invariants K_l are
                                  % undefined above l = 8.
CS_PHASE  = 0;                    % 0 == MRtrix's SH basis exactly. At SMI's
                                  % default of 1 the two differ by (-1)^m, a 180
                                  % degree rotation about z of every fODF.
PROTOCOL  = 'hcp_real_3shell.txt';
B0_SNAP   = 0.05;                 % ms/um^2. Any shell below this is treated as
                                  % exactly b = 0.
                                  %
                                  % The acquired .bval carries EFFECTIVE b, so
                                  % its "b = 0" volumes are b = 5 s/mm^2 rather
                                  % than 0 -- imaging gradients contribute a
                                  % little diffusion weighting of their own.
                                  % Left as acquired, that makes S(0)/S0 differ
                                  % slightly between kernels and puts a small
                                  % non-zero residual in Figure 1's lowest
                                  % shell, which is a distraction rather than a
                                  % result. Snapping it to 0 restores the exact
                                  % identity S(0)/S0 = 1 for every kernel.
                                  %
                                  % Only the b ~ 0 shell is snapped. The
                                  % diffusion-weighted shells keep their
                                  % acquired jitter (990-1005, 1985-2010,
                                  % 2980-3010 s/mm^2), so Step 1 still exercises
                                  % SMI's shell binning on real values.
                                  % Set to 0 to disable and use b as acquired.
NDIR_Q    = 3000;                 % quadrature directions for projecting a
                                  % sampled fODF onto plm
SEED      = 31415;                % RNG seed

% ---------------------------------------------- the constrained deconvolution
% Passed straight through to options.fODF_regularization. flag_nonneg = 1 is the
% arm being studied; it is OFF in the shipped toolbox defaults.
REG = struct('flag_nonneg', 1, 'lambda_tikhonov', 0.3);

% ------------------------------------------------------------ the CSD arms
% SSST-CSD and MSMT-CSD run inside this file, on the SAME |S_noisy| the SMI arm
% is fitted to, by calling the MRtrix3 binaries as subprocesses (Step 6b).
% Nothing here reimplements MRtrix: |dwi2fod|, |dwiextract| and |mrinfo| are the
% real thing, and the only MRtrix behaviour implemented locally is reading and
% writing its image format (|mrtrix_io.m|), which |mrinfo| checks on every run.
%
% *Why in this file rather than alongside it.* All three arms then see one
% simulation, one noise draw and one peak finder. There is no second
% implementation of the forward model to keep in step, and no question about
% whether the arms were compared on the same data -- they are the same array.
RUN_MRTRIX = 1;     % 0 scores the SMI arm alone, exactly as before

% *The response the CSD arms are given.* SMI has no response function, it has a
% kernel; the conversion is one line (Step 3). But which response is the honest
% one depends on what a "single fibre" is in this simulation, and here it is not
% a delta -- it is a Watson population at KAPPA.
%
%   'delta'      r_l(b) = K_l(b) sqrt((2l+1) 4 pi)
%                the exact analytic kernel response, sharper than real tissue.
%   'dispersed'  r_l(b) = K_l(b) p_l^Watson(KAPPA) sqrt((2l+1) 4 pi)
%                the kernel convolved with the same Watson the ground truth
%                disperses each population by -- the response a perfect
%                dwi2response run on this data would recover.
%
% Step 6b prints both next to the three real dwi2response estimators, so the
% choice is visible rather than buried. 'dispersed' is the default because a
% response estimated from real white matter has already absorbed dispersion, and
% giving CSD a sharper response than any estimator could produce would flatter
% it. Note the measured consequence, recorded in notebooks/README.md: MSMT-CSD's
% 60 degree error is response-limited and gets WORSE with the blunter response.
RESPONSE_MODE = 'dispersed';

% *WHICH TISSUE'S KERNEL THE CSD RESPONSE IS BUILT FROM. This is the assumption
% the manuscript is really about, so it is a knob and not a hard-coded choice.*
%
%   'healthy'  the healthy WM response is used for BOTH tissues (the default,
%              and what real CSD does)
%   'matched'  each tissue gets a response built from its own kernel
%
% *Why 'healthy' is the realistic setting.* A CSD response function is estimated
% once per subject -- or per study -- by selecting single-fibre white matter
% voxels and averaging them. Nobody estimates a separate response for edematous
% tissue: the edema is what is being imaged, not a population you can draw a
% reference from, and the selection heuristics (|dwi2response tournier|,
% |dhollander|) are built to find the most anisotropic voxels, which is exactly
% what edema is not. So on real data the edema voxels are deconvolved with a
% response averaged over HEALTHY white matter, and any mismatch between that
% response and the tissue actually present is an error the method has to live
% with.
%
% *Why that matters here more than anywhere else in this file.* It is the one
% place where the two arms are NOT on equal footing, and the inequality is real
% rather than an artefact:
%
% * SMI *estimates its kernel per voxel*. In edema it fits an edema kernel, so
%   it adapts. Step 6 already shows the estimate is imperfect, which is the
%   honest version of that advantage.
% * CSD *is given a fixed response*. In edema that response describes different
%   tissue from the one in the voxel.
%
% Setting this to 'matched' hands CSD a response it could not have on real data
% and makes the edema comparison flattering to it. That configuration is worth
% running as a CONTROL -- it separates "CSD is hurt by the response mismatch"
% from "CSD is hurt by the low anisotropic signal in edema" -- but it is not the
% configuration a manuscript claim should rest on.
%
% Step 6b prints the response it used against the response the tissue would have
% implied, so the size of the mismatch is on the record rather than assumed.
%
% *The full derivation of the response -- why it takes the form it does, what it
% assumes, and where each step is checked -- is in
% Reports/REPORT_CSD_response_derivation.md.* Read that before changing anything
% about how the response is built.
CSD_RESPONSE_KERNEL = 'healthy';

% *The two regularisation weights msmt_csd takes, and why they are set here
% rather than left at MRtrix's defaults.* This is the single most consequential
% setting in Step 6b and it was got wrong once, so it is a knob with its
% measurement attached rather than a silent default.
%
% |dwi2fod csd| and |dwi2fod msmt_csd| do NOT ship comparable defaults. The
% SSST algorithm's non-negativity constraint has strength 1; msmt_csd's
% |-neg_lambda| defaults to *1e-10*, i.e. essentially unregularised. Running
% both "at their defaults" therefore compares a constrained arm against an
% unconstrained one, and the MSMT fODF comes back much blunter -- l = 6 band
% power 0.31 against SSST's 1.11 -- which pulls the two lobes of a crossing
% together and displaces the peaks.
%
% Measured on the noise-free 60 degree crossing, peak separation (truth 60.00,
% band-limited truth at Lmax 6 gives 60.94):
%
%     true angle      SSST-CSD   MSMT default   MSMT neg_lambda 1
%        60 deg         60.94        48.67           60.94
%        75 deg         73.95        70.99           73.95
%        90 deg         89.36        86.38           89.36
%
% At MRtrix's default MSMT under-separates at EVERY angle, including 90 degrees,
% which is not a result any published MSMT-CSD comparison shows and was the tell
% that the setup was wrong rather than the method. |neg_lambda = 1| matches the
% SSST arm's constraint strength, which is the like-for-like choice.
%
% *This was NOT excluded by tuning to taste.* Ruled out first, each by direct
% measurement: the shells (msmt on b = 3 alone fails identically), the tissue
% count (WM-only identical), the response family (delta only partly helps), the
% peak finder (sh2peaks agrees to 0.06 deg), the SH basis (SSST is exact), and
% signal scaling (bit identical from S0 = 1 to 1e4).
%
% *Report these values with any MSMT number.* The sensitivity is large -- 48.67
% to 64.97 degrees across the range tested at a true 60 degree crossing -- so an
% MSMT result quoted without them is not reproducible.
MSMT_NEG_LAMBDA  = 1;       % MRtrix default 1e-10; 1 matches csd's constraint
MSMT_NORM_LAMBDA = 1e-3;    % MRtrix default 1e-10; see the note above

D_GM = 0.8;         % grey matter diffusivity, um^2/ms. Used ONLY to build the
                    % idealized isotropic response msmt_csd needs as its second
                    % tissue -- no simulated voxel contains grey matter, so what
                    % MSMT assigns there is a measurement of its own leakage.

% ------------------------------------------------------------------ assembled
% Collected into one struct so the rest of the file reads C.<knob> and there is
% exactly one place to change any of them.
C = struct();
C.D_FW     = D_FW;       C.KAPPA    = KAPPA;    % C.K_WM and C.PRESET are set
C.ANGLES   = ANGLES;     C.AXIS1    = AXIS1;    C.LMAX_GT  = LMAX_GT;
C.CS_PHASE = CS_PHASE;   C.PROTOCOL = PROTOCOL; C.NDIR_Q   = NDIR_Q;
C.SEED_MC  = SEED;                              % per kernel, in the loop
C.pick_grid      = MC.pick_grid;                 % shared utilities
C.rotate_about   = MC.rotate_about;
C.load_protocol  = @() MC.load_protocol_file(PROTOCOL);

% The fibre axes of each condition, built once from AXIS1 and the shared
% rotation so the geometry convention matches the pipeline exactly.
AX = cell(1, numel(ANGLES));
for ic = 1:numel(ANGLES)
    if ANGLES(ic) == 0
        AX{ic} = {AXIS1};
    else
        AX{ic} = {AXIS1, C.rotate_about(AXIS1, ANGLES(ic))};
    end
end
C.condition_axes = @(ic) AX{ic};

fprintf('\n=== SMI manuscript simulation ===\n');
if SMOKE_TEST
    fprintf('*** SMOKE_TEST = true: reduced sweep, INDICATIVE NUMBERS ONLY ***\n');
end
for ik = 1:NKERN
    fprintf('kernel %d/%d   : %-8s [f Da Depar Deperp fw] = %s, extra-axonal = %.2f\n', ...
            ik, NKERN, KERNELS{ik}.name, mat2str(KERNELS{ik}.K), ...
            1 - KERNELS{ik}.K(1) - KERNELS{ik}.K(5));
    fprintf('                %s\n', KERNELS{ik}.note);
end
fprintf('condition     : %g deg crossing only,  kappa = %g\n', ANGLES, KAPPA);
fprintf('sweep         : %d kernels x SNR %s x Lmax %s = %d fits of %d voxels each\n', ...
        NKERN, mat2str(SNR_LIST), mat2str(LMAX_LIST), ...
        NKERN*numel(SNR_LIST)*numel(LMAX_LIST), numel(ANGLES)*NREP);
fprintf('truth at Lmax %d (SMI''s kernel ceiling), CS_phase %d\n\n', LMAX_GT, CS_PHASE);

NSNR       = numel(SNR_LIST);
SIGMA_LIST = 1./SNR_LIST;                     % Inf -> 0
SNR_LABEL  = cell(1, NSNR);
for is = 1:NSNR
    if isinf(SNR_LIST(is))
        SNR_LABEL{is} = 'inf';
    else
        SNR_LABEL{is} = sprintf('%g', SNR_LIST(is));
    end
end
[~, SNR_ORD] = sort(SNR_LIST);  % ascending, Inf last: the order every figure
                                % and every plot uses, whatever order SNR_LIST
                                % was typed in

%% Which Lmax, and why the ground truth is stuck at 8
% This is the first thing to be clear about, because every result below depends
% on it and the numbers mean different things at different orders.
%
% *The fit runs at Lmax 4, 6 and 8, side by side.* Even spherical harmonics up
% to Lmax number |(Lmax/2+1)(Lmax+1)| = *15, 28 and 45* coefficients. The
% protocol has 90 directions per shell, so all three are comfortably determined
% by the data; the limit here is the model and the noise, not the sampling.
%
% *The ground truth is built at Lmax 8, and that is a ceiling, not a choice.*
% SMI's kernel rotational invariants |K_l(b)| are only defined up to l = 8
% (|SMI.RotInv_Kell_wFW_b_beta_TE_numerical|; ask for l = 10 and it errors).
% So the truth cannot be made sharper than the highest order being fitted.
%
% That has a consequence worth stating plainly:
%
% * At *Lmax 4 and 6* the truth carries angular detail the fit cannot represent,
%   which is the honest arrangement -- the fit is being asked to do something it
%   provably cannot do perfectly.
% * At *Lmax 8* the truth is exactly representable by the fit. The "ceiling"
%   reported in Step 7 is then no longer an independent bound, and any Lmax 8
%   result should be read as slightly flattering.
%
% The published comparison in |Reports/| ran every arm at *Lmax 6*, so 6 is the
% row to compare against the report. 4 is what the real-data driver |run_smi_batch_mod.m|
% uses. 8 is the report's own upper control.

fprintf('\n=== SMI manuscript simulation ===\n');
fprintf('protocol %s, NREP %d per condition per SNR, CS_phase %d\n', ...
        C.PROTOCOL, NREP, C.CS_PHASE);
fprintf('SNR sweep: %s   (%d noise levels)\n', strjoin(SNR_LABEL, ', '), NSNR);
fprintf('fitting at Lmax %s   (%s coefficients)\n', mat2str(LMAX_LIST), ...
        mat2str((LMAX_LIST/2+1).*(LMAX_LIST+1)));
fprintf('ground truth at Lmax %d, which is SMI''s kernel ceiling\n\n', LMAX_GT);

%% Step 1 -- the acquisition protocol
% A *real HCP 3-shell scheme*, supplied as FSL |.bval| / |.bvec| and tracked as
% text at |protocol/hcp_real_3shell.txt| so this walkthrough needs no Python.
% 288 volumes: 18 at b = 0 plus 90 each at nominal b = 1, 2, 3 ms/um^2.
%
% b is carried in ms/um^2 throughout. The |.bval| file is in scanner units, so
% the tracked table is the supplied values *divided by 1000*: b = 1 here is
% 1000 s/mm^2.
%
% Two properties of real acquisitions show up immediately, and neither is
% tidied away, because both are things code that only ever saw synthetic
% protocols would get wrong:
%
% * *The b = 0 volumes are not b = 0 as acquired* -- they are b = 5 s/mm^2,
%   because a .bval records EFFECTIVE b and the imaging gradients contribute a
%   little weighting of their own. |B0_SNAP| in the Configuration block sets
%   them to exactly 0, which restores the exact identity S(0)/S0 = 1 for every
%   kernel. They also carry unit direction vectors like every other volume,
%   even though the direction is meaningless there.
% * *The b values jitter within each shell* -- 18 distinct values across the
%   scheme. The forward model below uses the exact per-volume b; SMI bins them
%   into shells for the kernel fit, and Step 1 checks that it bins them the way
%   a human would.

% Read through the shared loader, which every other arm also uses, so there is
% one reader for the acquisition even though the filename is a knob above. It
% prints a loud warning if the .bvec is not unit -- expect one here, and see the
% note under Step 3 for why it matters.
[bvals, bvecs] = C.load_protocol();
Ndwi = numel(bvals);

% Treat the near-zero shell as exactly b = 0 (see B0_SNAP above).
n_snap = sum(bvals > 0 & bvals < B0_SNAP);
if n_snap > 0
    fprintf(['   note: %d volumes with 0 < b < %g snapped to exactly b = 0\n' ...
             '         (acquired as b = %g ms/um^2 = %.0f s/mm^2, effective b\n' ...
             '          from the imaging gradients)\n'], ...
            n_snap, B0_SNAP, max(bvals(bvals < B0_SNAP)), ...
            1000*max(bvals(bvals < B0_SNAP)));
    bvals(bvals < B0_SNAP) = 0;
end

fprintf('Step 1: %d volumes, %d distinct b values\n', Ndwi, numel(unique(bvals)));

% Let SMI group the shells, and report what it decided. This is the function the
% fit itself uses, and its shell assignment is used for every per-shell number
% below -- rather than a hand-rolled "within x% of the nominal b", which
% silently mixes neighbouring shells when they are a factor of ~1.5 apart.
[tbl, ~, shell_id] = SMI.Group_dwi_in_shells_b_beta_TE(bvals, [], [], []);
b_shell = tbl(1,:);
n_shell = tbl(3,:);
for i = 1:numel(b_shell)
    raw = bvals(shell_id == i);
    fprintf('   shell %d: b = %.3f ms/um^2, %3d directions, raw b in [%.3f, %.3f]\n', ...
            i, b_shell(i), n_shell(i), min(raw), max(raw));
end
dw = shell_id(:)' > 1;                        % everything above the b~0 shell

% CHECK 1. Every direction must be a unit vector. This passes trivially because
% the loader normalised them -- the number that matters is the one in the
% warning it printed just above, which is the error as the .bvec was supplied.
e_nrm = max(abs(sqrt(sum(bvecs.^2, 2)) - 1));
fprintf('   CHECK all directions are unit    max| |g|-1 | = %.2e   %s\n', ...
        e_nrm, VERDICT{1+(e_nrm < 1e-12)});

% CHECK 2. SMI must recover the four shells a human sees, from 18 distinct b
% values. If it split a shell, every rotational invariant downstream would be
% estimated from a fraction of the directions.
ok_shell = (numel(b_shell) == 4) && all(n_shell(:)' == [18 90 90 90]);
fprintf('   CHECK SMI bins %d b values into 4 shells [18 90 90 90]   %s\n', ...
        numel(unique(bvals)), VERDICT{1+ok_shell});

% CHECK 3. Within each diffusion-weighted shell the directions should cover the
% sphere evenly. A clumped scheme would make the l >= 4 coefficients unstable in
% a way that looks like a method failure.
fprintf('   nearest-neighbour direction spacing per shell:\n');
for i = 2:numel(b_shell)
    g = bvecs(shell_id == i, :);
    cc = abs(g*g'); cc(1:size(g,1)+1:end) = 0;
    nn = acosd(min(max(max(cc, [], 2), -1), 1));
    fprintf('     b = %.0f : mean %.2f deg, worst gap %.2f deg\n', ...
            b_shell(i), mean(nn), max(nn));
end
fprintf('\n');

%% Step 2 -- the ground truth fibre geometry
% Four conditions: a single fibre, and two equal fibre populations crossing at
% 30, 45 and 60 degrees. Each population is a *Watson* distribution with
% concentration |kappa = 16|, not a delta.
%
% That choice is deliberate and is the one most likely to be questioned, so:
% a response function estimated from real white matter has already absorbed
% fibre dispersion. A delta ground truth would create a response/truth mismatch
% that does not exist in practice, and would flatter whichever method sharpens
% most.
%
% The fODF is sampled on |NDIR_Q| quadrature directions and projected onto
% spherical harmonic coefficients |plm| in SMI's normalised convention
% |p_00 = 1|, in which *the fODF integrates to 1 in every voxel*.

dq    = H.dirs(C.NDIR_Q);
NCOND = numel(C.ANGLES);
COLHDR = cell(1, NCOND);
for ic = 1:NCOND
    if C.ANGLES(ic) == 0, COLHDR{ic} = 'single';
    else, COLHDR{ic} = sprintf('%d deg', C.ANGLES(ic)); end
end
L_gt  = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1)';

fodf_gt = zeros(size(dq,1), NCOND);
plm_gt  = zeros(NCOND, numel(L_gt)-1);
sh_gt   = zeros(NCOND, numel(L_gt));
axes_gt = cell(1, NCOND);
sep_deg = zeros(1, NCOND);

fprintf('Step 2: %d conditions\n', NCOND);
for ic = 1:NCOND
    axes_gt{ic} = C.condition_axes(ic);
    f = zeros(size(dq,1), 1);
    for k = 1:numel(axes_gt{ic})
        f = f + H.watson(dq, axes_gt{ic}{k}, C.KAPPA);
    end
    fodf_gt(:,ic) = f;
    p             = H.mixture_plm(f, dq, LMAX_GT, C.CS_PHASE);
    plm_gt(ic,:)  = p(:)';
    sh_gt(ic,:)   = ([1; p(:)] .* sqrt((2*L_gt+1)/(4*pi)))';

    if numel(axes_gt{ic}) == 2
        sep_deg(ic) = acosd(abs(axes_gt{ic}{1} * axes_gt{ic}{2}'));
    end
    fprintf('   condition %d: %2d deg, %d population(s), measured separation %.6f deg\n', ...
            ic, C.ANGLES(ic), numel(axes_gt{ic}), sep_deg(ic));
end

% CHECK 1. The angle between the two axes must be the angle we asked for. A
% wrong rotation axis would otherwise look plausible all the way to the end.
e_sep = max(abs(sep_deg - C.ANGLES));
fprintf('   CHECK crossing angles       max err  = %.2e deg   %s\n', ...
        e_sep, VERDICT{1+(e_sep < 1e-9)});

% CHECK 2. Reconstruct each fODF from its plm and integrate it over the sphere.
% In the p_00 = 1 convention the integral is 1 by construction, so anything
% else means the convention has been broken.
Yq     = SMI.get_even_SH(dq, LMAX_GT, C.CS_PHASE);
amp    = sh_gt * Yq';
mass   = mean(amp, 2) * 4*pi;                 % equal-area quadrature
e_mass = max(abs(mass - 1));
fprintf('   CHECK fODF integrates to 1  max|int-1| = %.2e   %s\n', ...
        e_mass, VERDICT{1+(e_mass < 1e-3)});

% CHECK 3. The fODF as *sampled* must be non-negative: it is a sum of Watson
% distributions, which are positive by construction. (These are unnormalised
% Watson amplitudes, so the floor is 1, not 0.)
mn_q = min(fodf_gt(:));
fprintf('   CHECK sampled fODF >= 0     min = %+.4f            %s\n', ...
        mn_q, VERDICT{1+(mn_q >= 0)});

% NOT a check, but the single most surprising number in this file, and it is
% better met here than discovered later:
%
% *the band-limited ground truth is itself negative in places.* Truncating a
% Watson mixture rings, exactly as truncating any Fourier series does, and the
% rings go below zero. So the "truth" the non-negativity constraint in Step 6 is
% asked to approach does not satisfy non-negativity either. This is why the
% constraint is a regularizer and not a statement of fact.
fprintf('   band-limited truth at Lmax %d -- peak, minimum, %% of sphere below 0:\n', LMAX_GT);
for ic = 1:NCOND
    fprintf('     %2d deg : peak %+.4f, min %+.4f, negative over %4.1f%% of directions\n', ...
            C.ANGLES(ic), max(amp(ic,:)), min(amp(ic,:)), 100*mean(amp(ic,:) < 0));
end
fprintf('   isotropic floor 1/(4*pi) = %.4f  (MRtrix iFOD2 default cutoff is 0.05)\n\n', ...
        1/(4*pi));

% ---- the Watson's own zonal coefficients, for the CSD arms' response
% These depend only on KAPPA and the band limit, not on the kernel, so they are
% computed once here rather than inside the kernel loop.
%
% A single fibre population in this simulation is a Watson at KAPPA, not a
% delta. Its normalised zonal coefficients p_l (with p_0 = 1) are what turn the
% delta response into the DISPERSION-MATCHED one in Step 6b:
%
%     r_l(b) = K_l(b) * p_l * sqrt((2l+1) * 4*pi)
%
% Computed by the same sampled-grid projection the ground truth uses, so the
% response and the truth cannot end up disagreeing about what a Watson is.
w1        = H.watson(dq, [0 0 1], C.KAPPA);
plm_w     = H.mixture_plm(w1, dq, LMAX_GT, C.CS_PHASE);
M_gt      = [];
for il = 0:2:LMAX_GT, M_gt = [M_gt, -il:il]; end %#ok<AGROW>
sh_w      = [1; plm_w(:)];
WATSON_PL = sh_w(M_gt(:) == 0);          % l = 0, 2, ..., LMAX_GT

% CHECK. The same numbers by 1-D quadrature, with no spherical harmonics
% anywhere: for an axially symmetric f(x), x = cos(theta),
%     p_l = int_0^1 f(x) P_l(x) dx / int_0^1 f(x) dx.
% Not identical -- the SH route is band limited at LMAX_GT and this one is not
% -- so the residual is the band-limiting error, not a mistake.
xq  = linspace(0, 1, 20001);
fq  = exp(C.KAPPA*xq.^2);
pl_q = zeros(size(WATSON_PL));
for il = 0:2:LMAX_GT
    switch il
        case 0, Pq = ones(size(xq));
        case 2, Pq = 1.5*xq.^2 - 0.5;
        case 4, Pq = (35*xq.^4 - 30*xq.^2 + 3)/8;
        case 6, Pq = (231*xq.^6 - 315*xq.^4 + 105*xq.^2 - 5)/16;
        case 8, Pq = (6435*xq.^8 - 12012*xq.^6 + 6930*xq.^4 - 1260*xq.^2 + 35)/128;
    end
    pl_q(il/2+1) = trapz(xq, fq.*Pq)/trapz(xq, fq);
end
e_pl = max(abs(WATSON_PL(:) - pl_q(:)));
fprintf('   Watson p_l at kappa = %g: %s\n', C.KAPPA, ...
        mat2str(round(WATSON_PL(:)'*1e5)/1e5));
fprintf('   CHECK Watson p_l, harmonics vs 1-D quadrature  max|err| = %.2e   %s\n', ...
        e_pl, VERDICT{1+(e_pl < 2e-4)});
fprintf('         (these blunt the CSD response in Step 6b; 1.0 would be a delta)\n\n');

%% Steps 3 to 8 run once per kernel
% Everything from here to the export is inside |for ik = 1:NKERN|, so the
% healthy and edema simulations run back to back off the same ground truth,
% the same protocol and the same noise seeds. Only the kernel differs.
%
% *These sections cannot be run as isolated Live Script cells* -- they are
% inside a loop, so run Steps 3 to 8 as a block. Each iteration stashes its
% results into |RUN{ik}| and the Figures section reads from there.

RUN = cell(1, NKERN);
for ik = 1:NKERN

C.K_WM   = KERNELS{ik}.K;
C.PRESET = KERNELS{ik}.name;
K        = C.K_WM;
fprintf('==================== kernel %d/%d: %s ====================\n', ...
        ik, NKERN, C.PRESET);

%% Step 3 -- the kernel, and the same kernel as a response function
% SMI has no response function. It has a *kernel*: the Standard Model
% compartment description |[f Da Depar Deperp fw]|, from which the rotational
% invariants |K_l(b)| follow analytically. CSD has the opposite arrangement, a
% non-parametric response estimated once and stored as zonal harmonics.
%
% They describe the same object. For a single fibre along z,
%
%  R(theta) = sum_l K_l(b) (2l+1) P_l(cos theta) = sum_l r_l Y_l0(theta)
%  r_l      = K_l(b) * sqrt((2l+1)*4*pi)
%
% and |r_l| is exactly one row of an MRtrix response |.txt| file. That identity
% is what lets an SMI kernel be handed to MRtrix at all, so it is checked here
% rather than assumed.

fprintf('Step 3: kernel [f Da Depar Deperp fw] = %s\n', mat2str(K));
fprintf('   zonal response r_l at the nominal shells -- this IS an MRtrix response file:\n');
fprintf('        b   ');
for l = 0:2:LMAX_GT, fprintf('%10s', sprintf('r_%d', l)); end
fprintf('\n');
b_nom = [0 1 2 3];
r_nom = RH.zh(K, b_nom, LMAX_GT, C.D_FW);
for i = 1:numel(b_nom)
    fprintf('     %4g   ', b_nom(i)); fprintf('%10.4f', r_nom(i,:)); fprintf('\n');
end
fprintf('   r_l falls fast with l: that decay is what makes deconvolution ill-conditioned,\n');
fprintf('   and it is why Lmax 8 is harder than Lmax 4 rather than simply better.\n');

% CHECK. Build the signal of a delta fODF along z two ways and compare.
%
% Route A, the zonal profile above. Route B, SMI's own forward model fed the
% plm of a delta. In the p_00 = 1 convention a delta along z is exact and needs
% no approximation: f_lm = Y_lm(z), which is zero for every m /= 0 and
% sqrt((2l+1)/4pi) at m = 0, so every p_l0 = 1 and every other p_lm = 0.
Mf = [];
for l = 2:2:LMAX_GT, Mf = [Mf, -l:l]; end
plm_delta = double(Mf(:)' == 0);

zax   = [0 0 1];
bsub  = bvals(dw);
gsub  = bvecs(dw,:);
S_fwd = H.signal(plm_delta, [K 1 1], bsub, ones(1,sum(dw)), zeros(1,sum(dw)), ...
                 gsub, LMAX_GT, C.CS_PHASE, C.D_FW);
theta = acos(max(min(gsub * zax', 1), -1));
r_each = RH.zh(K, bsub, LMAX_GT, C.D_FW);     % the exact per-volume b, jitter included
S_zon  = zeros(sum(dw), 1);
for i = 1:sum(dw)
    S_zon(i) = RH.profile(r_each(i,:), theta(i));
end
e_zon = max(abs(S_zon(:) - S_fwd(:)));
fprintf('   CHECK zonal response == SMI forward model  max|err| = %.2e   %s\n\n', ...
        e_zon, VERDICT{1+(e_zon < 1e-9)});

%% Step 4 -- forward convolution to a noise-free signal
% The signal is the fODF convolved with the kernel. In spherical harmonics
% convolution is a product:
%
%  S(u)/S0 = sum_lm K_l(b) p_lm Y_lm(u) sqrt((2l+1)*4*pi)
%
% This is exactly the expression |SMI.get_plm_from_S_and_kernel| inverts. Using
% it here is less circular than it looks -- the test is whether the *inverse*
% recovers the input, which is Step 6 -- but to make sure the harmonic
% machinery itself is right, the same signal is also computed with no spherical
% harmonics anywhere, by direct numerical convolution over the quadrature grid:
%
%  S(u) = sum_q w_q * fODF(n_q) * K(u.n_q),   with sum_q w_q = 1

S_clean = zeros(NCOND, Ndwi);
for ic = 1:NCOND
    s = H.signal(plm_gt(ic,:), [K 1 1], bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                 bvecs, LMAX_GT, C.CS_PHASE, C.D_FW);
    S_clean(ic,:) = s(:)';
end

f_ = K(1); Da_ = K(2); Dep_ = K(3); Dpp_ = K(4); fw_ = K(5);
cosang  = bvecs * dq';
Kmat    = f_         * exp(-bvals' .* (Da_*cosang.^2)) + ...
          (1-f_-fw_) * exp(-bvals' .* (Dpp_ + (Dep_-Dpp_)*cosang.^2)) + ...
          fw_        * exp(-bvals' * C.D_FW);
S_exact = zeros(NCOND, Ndwi);
for ic = 1:NCOND
    w = fodf_gt(:,ic) / sum(fodf_gt(:,ic));
    S_exact(ic,:) = (Kmat * w)';
end

fprintf('Step 4: noise-free signal\n');
fprintf('   mean signal per shell (shells as SMI assigned them in Step 1):\n');
fprintf('        b   ');  fprintf('%10s', COLHDR{:}); fprintf('\n');
for i = 1:numel(b_shell)
    fprintf('     %4.2f   ', b_shell(i));
    fprintf('%10.4f', mean(S_clean(:, shell_id == i), 2)); fprintf('\n');
end

% CHECK 1. The b~0 signal must be essentially 1, because sigma is set to 1/SNR
% below and that only means the requested SNR if S0 = 1. It is not exactly 1
% here: this protocol's "b = 0" is b = 5 s/mm^2, which attenuates by a few parts
% per thousand. That is a property of the scan, so it is measured, not assumed.
S_b0 = mean(S_clean(:, ~dw), 2);
e_s0 = max(abs(S_b0 - 1));
if B0_SNAP > 0
    fprintf('   CHECK S(b=0) == 1 exactly        max|S-1| = %.2e   %s\n', ...
            e_s0, VERDICT{1+(e_s0 < 1e-12)});
else
    fprintf('   CHECK S(b~0) is within 1%% of 1   max|S-1| = %.2e   %s\n', ...
            e_s0, VERDICT{1+(e_s0 < 0.01)});
end

% CHECK 2. Harmonics against direct convolution. These do NOT agree to machine
% precision, and should not: the harmonic route is band limited at LMAX_GT and
% the direct sum is not. The residual below IS the band-limiting error of the
% ground truth, and Step 5 puts it next to the noise.
e_sh = max(abs(S_clean(:) - S_exact(:)));
fprintf('   CHECK harmonics vs direct convolution  max|err| = %.2e\n', e_sh);
fprintf('         (band limiting, not an error -- it falls as LMAX_GT rises)\n\n');

%% Step 5 -- Rician noise, one block of voxels per SNR
% Complex Gaussian noise is added to a real signal and the magnitude taken,
% which is exactly Rician:
%
%  S_noisy = sqrt( (S + sigma*n1)^2 + (sigma*n2)^2 ),   n1, n2 ~ N(0,1)
%
% |sigma = 1/SNR|, with SNR defined against the b = 0 signal.
%
% *The sweep is laid out as one contiguous block of voxels per SNR*, each block
% holding all |NCOND| conditions at |NREP| realisations, in the same voxel order
% |gen_montecarlo.m| uses. Step 6 then hands one block at a time to |SMI.fit|,
% so every fit sees exactly one noise level -- which is what the pipeline does,
% and it keeps each fit's memory to one block rather than the whole sweep.
%
% *Each SNR is seeded separately*, from |SEED| in the Configuration block offset by the SNR
% index. Two reasons: each noise level is reproducible on its own, and adding or
% removing an entry from |SNR_LIST| does not silently change the realisations of
% every other entry. The cost is that the blocks are independent rather than
% sharing common random numbers, so a difference between two SNRs carries the
% Monte Carlo error of both.

NVOX_SNR = NCOND*NREP;                  % voxels in one SNR block: what a fit sees
NVOX     = NVOX_SNR*NSNR;               % voxels in the whole sweep
GRID_SNR = C.pick_grid(NVOX_SNR);       % the 3D grid one SMI.fit call is given
cond_id  = repmat(repelem((1:NCOND)', NREP, 1), NSNR, 1);
snr_id   = repelem((1:NSNR)', NVOX_SNR, 1);
S_rep    = S_clean(cond_id, :);
S_noisy  = zeros(NVOX, Ndwi);

for is = 1:NSNR
    rows = find(snr_id == is);
    sg   = SIGMA_LIST(is);
    rand('seed', C.SEED_MC + is); randn('seed', C.SEED_MC + is);
    Sb   = S_rep(rows, :);
    % At sigma = 0 this is sqrt(S^2) = S, i.e. the noise-free signal exactly.
    % The randn draws still happen so the noise-free arm costs the same and the
    % stream is the same shape at every SNR.
    S_noisy(rows,:) = sqrt((Sb + sg*randn(size(Sb))).^2 + ...
                           (     sg*randn(size(Sb))).^2);
end

fprintf('Step 5: %d conditions x %d reps x %d SNR = %d voxels\n', ...
        NCOND, NREP, NSNR, NVOX);
fprintf('   each SMI.fit sees one SNR block: %d voxels, grid %s\n', ...
        NVOX_SNR, mat2str(GRID_SNR));

% CHECK 1. Recover sigma from the simulated data at every SNR, rather than
% trusting the value typed in. This is an independent read of the noise level,
% and at SNR = inf it is the statement that the "noisy" signal is bit identical
% to the noise-free one.
%
% *Read at b = 0, not over all volumes, and the difference is not cosmetic.*
% Rician noise is a magnitude, so the residual S_noisy - S_clean has standard
% deviation sigma only where the signal is well ABOVE the noise floor. It is
% not, in the diffusion-weighted shells of a low-signal tissue: the edema kernel
% leaves a mean b = 3 signal of 0.0399 against sigma = 0.1, i.e. the signal sits
% two and a half times BELOW the noise, and there the magnitude operation
% compresses the residual and this estimate reads ~14% low. That is the Rician
% floor -- the same effect the table further down measures directly -- and not
% an error in the simulation.
%
% At b = 0 the signal is exactly 1 by construction, so S/sigma is exactly the
% requested SNR and this is the one place in the data where the residual is
% clean Gaussian. Both numbers are printed, because the gap between them is
% itself the Rician floor and is worth seeing.
b0_vol = ~dw;
for is = 1:NSNR
    rows  = find(snr_id == is);
    rblk  = S_noisy(rows,:) - S_rep(rows,:);
    rb0   = rblk(:, b0_vol);
    s_hat = std(rb0(:));                 % at b = 0, where S = 1 exactly
    s_all = std(rblk(:));                % over every volume, for comparison
    sg    = SIGMA_LIST(is);
    if sg == 0
        ok_is = (max(abs(rblk(:))) == 0);
        fprintf('   CHECK SNR %-4s  noise free, residual is exactly zero        %s\n', ...
                SNR_LABEL{is}, VERDICT{1+ok_is});
    else
        rel   = abs(s_hat - sg)/sg;
        ok_is = (rel < 0.10);
        fprintf(['   CHECK SNR %-4s  recovered sigma %.5f vs %.5f, %4.1f%% off  %s' ...
                 '   (all volumes: %.5f)\n'], ...
                SNR_LABEL{is}, s_hat, sg, 100*rel, VERDICT{1+ok_is}, s_all);
    end
end

% Where the band limit sits relative to the noise. At a finite SNR the noise is
% far larger than the truth's band-limiting error, so the band limit is not what
% limits the result. At SNR = inf it is the ONLY error left, which is exactly
% what the ceiling line of Step 7 measures.
sg_min = min(SIGMA_LIST(SIGMA_LIST > 0));
if isempty(sg_min)
    fprintf('   every SNR in the sweep is noise free, so the truth''s band-limiting\n');
    fprintf('   error %.2e is the only error anywhere below\n', e_sh);
else
    fprintf('   smallest non-zero sigma %.4f is %.0fx the band-limiting error %.2e,\n', ...
            sg_min, sg_min/e_sh, e_sh);
    fprintf('   so noise dominates the band limit at every finite SNR in the sweep\n');
end

% CHECK 2. A magnitude is strictly positive, which is why SMI is told
% NoiseBias = 'Rician' in Step 6.
mn = min(S_noisy(:));
fprintf('   CHECK signal strictly positive   min = %.4f            %s\n', ...
        mn, VERDICT{1+(mn > 0)});

% The Rician floor, per SNR. This is the bias the magnitude operation puts into
% the highest shell, and watching it grow as SNR falls is the reason
% NoiseBias = 'Rician' is not optional at the bottom of the sweep.
hib = (shell_id(:)' == numel(b_shell));
fprintf('   Rician floor at b = %.0f, mean over the shell:\n', b_shell(end));
fprintf('        SNR     noisy   noise free      upward bias\n');
for is = 1:NSNR
    rows = find(snr_id == is);
    mnz  = mean(mean(S_noisy(rows,hib)));
    mcl  = mean(mean(S_rep(rows,hib)));
    fprintf('     %6s  %8.4f     %8.4f     %+8.4f\n', SNR_LABEL{is}, mnz, mcl, mnz-mcl);
end

% NOT a check, but the thing about a sweep that is easiest to get wrong:
% *SMI does not fit the kernel at the sigma you pass in.* It normalises sigma by
% the measured b = 0 signal, bins the result into Nlevels equal bins spanning
% sigma_norm_limits, and trains one polynomial regression per occupied bin,
% evaluated at the bin CENTRE (SMI.m:2222-2302). Two consequences only a sweep
% makes visible: two SNRs that land in the same bin share a regression trained
% at one noise level, and an SNR whose sigma exceeds the top of the range is
% clamped and fitted with a regression trained at a LOWER noise level than the
% data actually has. Both constants are hard coded in SMI.m, so they are
% restated here rather than read back, and this block needs updating if they
% change.
SMI_SIGMA_LIMITS = [0 0.2];             % SMI.m:562
SMI_NLEVELS      = 10;                  % SMI.m:388, the default Nlevels
lev_w = (SMI_SIGMA_LIMITS(2)-SMI_SIGMA_LIMITS(1))/SMI_NLEVELS;
fprintf('   how SMI will bin these noise levels (Nlevels %d over sigma_norm %s):\n', ...
        SMI_NLEVELS, mat2str(SMI_SIGMA_LIMITS));
fprintf('        SNR     sigma   bin   trained at\n');
for is = 1:NSNR
    lev = min(max(floor(SIGMA_LIST(is)/lev_w) + 1, 1), SMI_NLEVELS);
    fprintf('     %6s  %8.4f  %4d     %8.4f', SNR_LABEL{is}, SIGMA_LIST(is), ...
            lev, (lev-0.5)*lev_w);
    if SIGMA_LIST(is) > SMI_SIGMA_LIMITS(2)
        fprintf('   <-- ABOVE the trained range, clamped');
    end
    fprintf('\n');
end
fprintf('   These bins are nominal: sigma is normalised by the MEASURED b = 0\n');
fprintf('   signal, so voxels of one SNR can straddle a bin edge.\n\n');

%% Step 6 -- SMI.fit, at each Lmax in turn
% The real toolbox at its shipped defaults, with the constrained deconvolution
% turned on. Two options are not defaults and both matter:
%
% * |CS_phase = 0|. At SMI's default of 1 the spherical harmonic basis differs
%   from MRtrix's by |(-1)^m|, which for even l is a 180 degree rotation about
%   z of every fODF -- 71.5 degrees of peak error, verified against MRtrix's
%   own |sh2peaks| in |check_mrtrix_basis.sh|. At 0 the two bases are identical,
%   which is what makes the Step 8 export directly loadable.
% * |fODF_regularization.flag_nonneg = 1|. Off by default in the toolbox; this
%   is the arm being studied. |lambda_nonneg| is left at its shipped value and
%   printed back as confirmation.
%
% The whole fit is repeated for each Lmax *and each SNR*, so everything
% downstream can be read as a function of angular order and of noise. Each call
% is given one SNR block and the scalar sigma of that block, which is what
% |gen_montecarlo.m| does -- no fit ever sees two noise levels at once.

fits = cell(1, numel(LMAX_LIST));
fprintf('Step 6: %d SMI.fit calls (%d Lmax x %d SNR), %d voxels each\n', ...
        numel(LMAX_LIST)*NSNR, numel(LMAX_LIST), NSNR, NVOX_SNR);
for iL = 1:numel(LMAX_LIST)
    Lf    = LMAX_LIST(iL);
    Lv    = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
    ncoef = numel(Lv);

    sh_all    = zeros(NVOX, ncoef);
    kern_all  = [];
    conv_all  = false(NVOX, 1);
    nonfinite = 0;
    sec_L     = 0;

    for is = 1:NSNR
        rows  = find(snr_id == is);
        dwi_b = reshape(S_noisy(rows,:), [GRID_SNR Ndwi]);

        options = struct();
        options.b     = bvals;
        options.dirs  = bvecs;
        options.sigma = SIGMA_LIST(is)*ones(GRID_SNR);
        options.mask  = true(GRID_SNR);
        options.compartments  = {'IAS','EAS','FW'};
        options.NoiseBias     = 'Rician';
        options.Lmax          = [0 Lf Lf Lf];
        options.CS_phase      = C.CS_PHASE;
        options.D_FW          = C.D_FW;
        options.flag_fit_fODF = 1;
        options.fODF_regularization = REG;   % from the Configuration block

        t0  = tic;
        out = SMI.fit(dwi_b, options);
        el  = toc(t0);
        sec_L = sec_L + el;

        plm = reshape(out.plm, [NVOX_SNR ncoef-1]);
        sh  = [ones(NVOX_SNR,1) plm] .* repmat(sqrt((2*Lv+1)/(4*pi))', NVOX_SNR, 1);
        nonfinite = nonfinite + sum(~isfinite(sh(:)));
        sh(~isfinite(sh)) = 0;
        sh_all(rows,:) = sh;

        kb = reshape(out.kernel, [NVOX_SNR size(out.kernel,4)]);
        if isempty(kern_all)
            kern_all = zeros(NVOX, size(kb,2));
        end
        kern_all(rows,:) = kb;
        conv_all(rows)   = (out.fODF_regularization.flag_converged(:) == 1);

        fprintf('   Lmax %d, SNR %-4s: %2d coefficients, %6.1f s  (lambda_nonneg = %g, lambda_tikhonov = %g)\n', ...
                Lf, SNR_LABEL{is}, ncoef, el, ...
                out.fODF_regularization.lambda_nonneg, ...
                out.fODF_regularization.lambda_tikhonov);
    end

    fits{iL} = struct('Lmax', Lf, 'sh', sh_all, 'nonfinite', nonfinite, ...
                      'kernel', kern_all, 'converged', conv_all, ...
                      'seconds', sec_L);
end

% The kernel is estimated from rotational invariants, which do not depend on
% the fODF's Lmax in the same way -- so this table should be nearly flat across
% Lmax. If it is not, the kernel fit is being disturbed by the fODF fit. Down
% the SNR axis it is NOT expected to be flat: this is where the noise shows up
% first, well before the fODF does.
knm = {'f','Da','Depar','Deperp','fw'};
fprintf('   kernel recovery, median over the %d voxels of each block (truth in brackets):\n', ...
        NVOX_SNR);
fprintf('        Lmax   SNR  ');
for j = 1:5, fprintf('%16s', sprintf('%s [%.2f]', knm{j}, K(j))); end
fprintf('\n');
for iL = 1:numel(LMAX_LIST)
    for k = 1:NSNR
        is   = SNR_ORD(k);
        rows = find(snr_id == is);
        fprintf('        %4d  %4s  ', LMAX_LIST(iL), SNR_LABEL{is});
        for j = 1:5
            col = fits{iL}.kernel(rows,j); col = col(isfinite(col));
            fprintf('%16.3f', median(col));
        end
        fprintf('\n');
    end
end
fprintf(['   The kernel is NOT recovered exactly, and the fODF below is deconvolved with\n' ...
         '   this ESTIMATED kernel rather than the true one. That is deliberate: it is what\n' ...
         '   happens on real data. A 3-shell linear-encoding protocol constrains the kernel\n' ...
         '   only loosely -- see README.md "Useful tips".\n']);

% CHECK. The constrained deconvolution must converge in every voxel at every
% Lmax, and no NaN may reach the fODF: a NaN in an SH volume breaks downstream
% tractography silently.
all_conv = true; all_fin = true;
for iL = 1:numel(LMAX_LIST)
    all_conv = all_conv && all(fits{iL}.converged == 1);
    all_fin  = all_fin  && (fits{iL}.nonfinite == 0);
end
fprintf('   CHECK deconvolution converged everywhere               %s\n', VERDICT{1+all_conv});
% Per SNR, because if convergence is going to fail it will fail at the bottom of
% the sweep, and an aggregate "ok" over 12000 voxels would hide which block.
fprintf('   converged fraction per block:\n');
fprintf('        Lmax  ');
for k = 1:NSNR, fprintf('%8s', SNR_LABEL{SNR_ORD(k)}); end
fprintf('\n');
for iL = 1:numel(LMAX_LIST)
    fprintf('        %4d  ', LMAX_LIST(iL));
    for k = 1:NSNR
        rows = find(snr_id == SNR_ORD(k));
        fprintf('%7.1f%%', 100*mean(fits{iL}.converged(rows)));
    end
    fprintf('\n');
end
fprintf('   CHECK fODF all finite at every Lmax                    %s\n', VERDICT{1+all_fin});

% CHECK. The p_00 = 1 convention must survive the fit at every Lmax.
e_l0 = 0;
for iL = 1:numel(LMAX_LIST)
    e_l0 = max(e_l0, max(abs(fits{iL}.sh(:,1) - 1/sqrt(4*pi))));
end
fprintf('   CHECK p_00 == 1 convention held   max|err| = %.2e   %s\n\n', ...
        e_l0, VERDICT{1+(e_l0 < 1e-12)});

%% Step 6b -- the CSD arms: dwi2fod csd and dwi2fod msmt_csd
% *The same noisy signal Step 5 built*, written out as an MRtrix image and
% deconvolved by the real MRtrix3 binaries. Not a reimplementation, and not a
% separate simulation: |S_noisy| is the array |SMI.fit| was just given, so the
% three arms differ only in the deconvolution.
%
% This is what makes the comparison answerable. There is no second forward
% model to keep in step, no second noise draw, and no cross-language
% reconciliation -- the arms share the same voxels by construction.
%
% Four things have to be right, and each gets a CHECK rather than an assumption:
%
% # *MRtrix must bin the shells the way SMI did.* The acquired b values jitter
%   (18 distinct values across 4 shells) and MRtrix does its own clustering. The
%   shells are read back with |mrinfo -shell_bvalues|, compared against SMI's
%   binning, and *the response is evaluated at MRtrix's own average b per
%   shell* -- not at the nominal 0/1/2/3, which would attach each response row
%   to a b nobody acquired at.
% # *The response has to survive the round trip.* It is written, read back off
%   disk, and compared against the array that was written.
% # *Lmax has to match, and the response file is what enforces it.* |-lmax| is
%   NOT passed. MRtrix's documented default is "the lmax of the corresponding
%   response function, based on its number of coefficients, up to a maximum
%   of 8" -- and one response file is written per Lmax with exactly |Lf/2+1|
%   columns, so the order is already pinned by the file. Verified rather than
%   assumed: with and without the flag, |dwi2fod| returns the same 15, 28 and
%   45 coefficients at Lmax 4, 6 and 8, and the FODs are bit identical. The
%   same mechanism gives msmt_csd's GM and CSF lmax 0, since their responses
%   are one column wide.
% # *Scale is not comparable, and does not need to be.* An SMI fODF has
%   |p_00 = 1| and integrates to 1; an MRtrix FOD is unnormalised and its
%   amplitude carries apparent fibre density. Step 7's peak finder subtracts
%   *each voxel's own l = 0 term* rather than the constant |1/(4*pi)|, which
%   makes it scale free. For an SMI fODF those are the same number, and that is
%   checked.
%
% *MSMT-CSD needs three tissues and this simulation contains one.* Its second
% and third responses are idealized isotropic ones at |D_GM| and |C.D_FW|. No
% simulated voxel contains grey matter or free-standing CSF, so what MSMT
% assigns to those compartments is a measurement of its own leakage, printed
% below rather than discarded.

% Every arm carries the same two things: a name and one [NVOX x ncoef] SH matrix
% per Lmax. Step 7 does not know or care which arm produced which matrix.
ARMS = {};
ARMS{end+1} = struct('name', 'SMI', 'sh', {cellfun(@(f) f.sh, fits, ...
                                                   'UniformOutput', false)});

if ~RUN_MRTRIX
    fprintf('Step 6b: SKIPPED (RUN_MRTRIX = 0), scoring the SMI arm alone\n\n');
else
MR   = mrtrix_io();
MDIR = fullfile(pkgdir, 'mrtrix');
if ~exist(MDIR, 'dir'), mkdir(MDIR); end
TAG  = sprintf('ms_%s', C.PRESET);
fprintf('Step 6b: the CSD arms, via MRtrix3, in %s\n', MDIR);

% ---- write the DWI the binaries read
% b in s/mm^2, which is what MRtrix expects; this file carries ms/um^2
% everywhere else. Voxel order is column major over GRID_ALL, so voxel v of
% S_noisy is voxel v of the image and Step 7 indexes them identically.
GRID_ALL = C.pick_grid(NVOX);
f_dwi    = fullfile(MDIR, [TAG '_dwi']);
f_mask   = fullfile(MDIR, [TAG '_mask']);
MR.write([f_dwi '.mif'],  reshape(S_noisy, [GRID_ALL Ndwi]), ...
         struct('grad', [bvecs bvals(:)*1000]));
MR.write([f_mask '.mif'], ones(GRID_ALL), struct('datatype', 'UInt8'));
fprintf('   wrote %s_dwi.mif  [%s x %d] with dw_scheme\n', ...
        TAG, mat2str(GRID_ALL), Ndwi);

% CHECK. mrinfo must agree with what we think we wrote. This is the one place
% the locally implemented image writer is checked against MRtrix itself.
[st, txt] = system(sprintf('mrinfo -size "%s.mif" 2>&1', f_dwi));
sz_mr = sscanf(txt, '%f')';
ok_sz = (st == 0) && numel(sz_mr) == 4 && all(sz_mr == [GRID_ALL Ndwi]);
fprintf('   CHECK mrinfo reads back size %s               %s\n', ...
        mat2str(sz_mr), VERDICT{1+ok_sz});
if ~ok_sz
    fprintf(2, '%s\n', txt);
    error(['MRtrix could not read the image this file wrote. Is mrtrix3 on ' ...
           'the PATH? Set RUN_MRTRIX = 0 to score the SMI arm alone.']);
end

% ---- what MRtrix thinks the shells are
[~, txt] = system(sprintf('mrinfo -shell_bvalues "%s.mif" 2>/dev/null', f_dwi));
b_mr = sscanf(txt, '%f')';
[~, txt] = system(sprintf('mrinfo -shell_sizes "%s.mif" 2>/dev/null', f_dwi));
n_mr = sscanf(txt, '%f')';
ok_shells = (numel(b_mr) == numel(b_shell)) && all(n_mr(:)' == n_shell(:)') && ...
            all(diff(b_mr) > 0);
fprintf('   CHECK MRtrix and SMI agree on the shells (%d, sizes %s)   %s\n', ...
        numel(b_mr), mat2str(n_mr), VERDICT{1+ok_shells});
fprintf('        shell     SMI b        MRtrix b     difference\n');
for i = 1:numel(b_mr)
    fprintf('        %4d   %9.5f   %13.5f   %11.2e\n', ...
            i, b_shell(i), b_mr(i)/1000, b_shell(i) - b_mr(i)/1000);
end

% ---- the single-shell subset for SSST-CSD
% dwi2fod csd is single shell and is run on the top shell, which is where SSST
% CSD is normally run. dwiextract is given MRtrix's own shell b value rather
% than the nominal 3000, so it cannot miss a jittered volume.
f_b3 = fullfile(MDIR, [TAG '_b3']);
[st, txt] = system(sprintf( ...
    'dwiextract "%s.mif" -shells 0,%g "%s.mif" -force -quiet 2>&1', ...
    f_dwi, b_mr(end), f_b3));
if st ~= 0, fprintf(2, '%s\n', txt); error('dwiextract failed'); end

% ---- the responses, one set per Lmax
% *Which kernel the response is built from is a setting, not this iteration's
% kernel.* See CSD_RESPONSE_KERNEL in the Configuration block. At the default,
% 'healthy', both tissues are deconvolved with the healthy WM response, because
% that is what a response estimated by population-averaging single-fibre white
% matter actually is. The edema arm is therefore deconvolved with a response
% describing DIFFERENT TISSUE, which is the situation on real data and is the
% asymmetry the manuscript is about: SMI re-estimates its kernel per voxel and
% adapts, CSD cannot.
if strcmp(CSD_RESPONSE_KERNEL, 'matched')
    K_RESP    = K;
    resp_name = C.PRESET;
else
    knames = cellfun(@(s) s.name, KERNELS, 'UniformOutput', false);
    ir = find(strcmp(knames, CSD_RESPONSE_KERNEL), 1);
    if isempty(ir)
        error(['CSD_RESPONSE_KERNEL = ''%s'' names no kernel. ' ...
               'Available: %s, or ''matched''.'], ...
              CSD_RESPONSE_KERNEL, strjoin(knames, ', '));
    end
    K_RESP    = KERNELS{ir}.K;
    resp_name = KERNELS{ir}.name;
end

r_delta_top = RH.zh(K_RESP, b_mr(end)/1000, LMAX_GT, C.D_FW);
r_disp_top  = r_delta_top .* WATSON_PL(:)';
fprintf('   the response, normalised to l = 0, at the top shell b = %.2f:\n', ...
        b_mr(end)/1000);
fprintf('        %-22s', 'delta (exact kernel)');
fprintf('%9.4f', r_delta_top/r_delta_top(1)); fprintf('\n');
fprintf('        %-22s', 'dispersion matched');
fprintf('%9.4f', r_disp_top/r_disp_top(1)); fprintf('\n');
fprintf('        RESPONSE_MODE = ''%s'',  built from the ''%s'' kernel %s\n', ...
        RESPONSE_MODE, resp_name, mat2str(K_RESP));

% NOT a check, a measurement, and one of the numbers the manuscript wants: how
% far the response the CSD arms were GIVEN is from the response this tissue
% would actually imply. Zero for the healthy arm at the default; the edema arm
% is where it bites, and its size is the response mismatch CSD carries into
% edema on real data.
r_tissue_top = RH.zh(K, b_mr(end)/1000, LMAX_GT, C.D_FW);
if strcmp(RESPONSE_MODE, 'dispersed')
    r_tissue_top = r_tissue_top .* WATSON_PL(:)';
end
n_given  = r_disp_top/r_disp_top(1);
n_tissue = r_tissue_top/r_tissue_top(1);
if max(abs(n_given - n_tissue)) < 1e-12
    fprintf('        response matches this tissue exactly (no mismatch)\n');
else
    fprintf('        MISMATCH -- this tissue (%s) would imply', C.PRESET);
    fprintf('%9.4f', n_tissue); fprintf('\n');
    fprintf('        difference                              ');
    fprintf('%9.4f', n_given - n_tissue); fprintf('\n');
    fprintf(['        The %s arm is deconvolved with the %s response. That is\n' ...
             '        deliberate and it is what real CSD does -- a response is\n' ...
             '        estimated from single-fibre WM, never from edema. Set\n' ...
             '        CSD_RESPONSE_KERNEL = ''matched'' for the control.\n'], ...
            C.PRESET, resp_name);
end

f_resp = cell(1, numel(LMAX_LIST));
for iL = 1:numel(LMAX_LIST)
    Lf   = LMAX_LIST(iL);
    r_wm = RH.zh(K_RESP, b_mr/1000, Lf, C.D_FW);       % [nshell x (Lf/2+1)]
    if strcmp(RESPONSE_MODE, 'dispersed')
        r_wm = r_wm .* repmat(WATSON_PL(1:Lf/2+1)', size(r_wm,1), 1);
    end
    r_gm  = exp(-(b_mr(:)/1000)*D_GM)*sqrt(4*pi);      % [nshell x 1]
    r_csf = exp(-(b_mr(:)/1000)*C.D_FW)*sqrt(4*pi);

    R = struct();
    R.wm    = fullfile(MDIR, sprintf('%s_resp_wm_lmax%d.txt', TAG, Lf));
    R.gm    = fullfile(MDIR, sprintf('%s_resp_gm.txt',  TAG));
    R.csf   = fullfile(MDIR, sprintf('%s_resp_csf.txt', TAG));
    R.wm_b3 = fullfile(MDIR, sprintf('%s_resp_wm_b3_lmax%d.txt', TAG, Lf));
    RH.write_response(R.wm,    r_wm);
    RH.write_response(R.gm,    r_gm);
    RH.write_response(R.csf,   r_csf);
    RH.write_response(R.wm_b3, r_wm(end,:));           % single shell, one row
    f_resp{iL} = R;

    % CHECK. The file MRtrix will read must hold the numbers computed here.
    % A disk round trip, not a comparison of an array with itself.
    e_rt = max(max(abs(RH.read_response(R.wm) - r_wm)));
    fprintf('   CHECK Lmax %d response round trip  max|err| = %.2e   %s\n', ...
            Lf, e_rt, VERDICT{1+(e_rt < 1e-7)});
end

% ---- run the binaries
SH_CSD  = cell(1, numel(LMAX_LIST));
SH_MSMT = cell(1, numel(LMAX_LIST));
tissue_share = nan(numel(LMAX_LIST), 3);

for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);
    R  = f_resp{iL};
    nc = (Lf/2+1)*(Lf+1);

    f_msmt = fullfile(MDIR, sprintf('%s_msmtfod_lmax%d', TAG, Lf));
    f_mgm  = fullfile(MDIR, sprintf('%s_msmtgm_lmax%d',  TAG, Lf));
    f_mcsf = fullfile(MDIR, sprintf('%s_msmtcsf_lmax%d', TAG, Lf));
    f_csd  = fullfile(MDIR, sprintf('%s_csdfod_lmax%d',  TAG, Lf));

    t0 = tic;
    [st, txt] = system(sprintf( ...
        ['dwi2fod msmt_csd "%s.mif" "%s" "%s.mif" "%s" "%s.mif" "%s" "%s.mif" ' ...
         '-mask "%s.mif" -neg_lambda %g -norm_lambda %g -force -quiet 2>&1'], ...
        f_dwi, R.wm, f_msmt, R.gm, f_mgm, R.csf, f_mcsf, f_mask, ...
        MSMT_NEG_LAMBDA, MSMT_NORM_LAMBDA));
    if st ~= 0, fprintf(2, '%s\n', txt); error('dwi2fod msmt_csd failed at Lmax %d', Lf); end
    t_msmt = toc(t0);

    t0 = tic;
    [st, txt] = system(sprintf( ...
        ['dwi2fod csd "%s.mif" "%s" "%s.mif" -mask "%s.mif" ' ...
         '-force -quiet 2>&1'], f_b3, R.wm_b3, f_csd, f_mask));
    if st ~= 0, fprintf(2, '%s\n', txt); error('dwi2fod csd failed at Lmax %d', Lf); end
    t_csd = toc(t0);

    Vm = MR.read([f_msmt '.mif']); SH_MSMT{iL} = reshape(Vm, [NVOX size(Vm,4)]);
    Vc = MR.read([f_csd  '.mif']); SH_CSD{iL}  = reshape(Vc, [NVOX size(Vc,4)]);
    Vg = MR.read([f_mgm  '.mif']); Vf = MR.read([f_mcsf '.mif']);

    fprintf('   Lmax %d: msmt_csd %.1f s, csd %.1f s, %d coefficients each\n', ...
            Lf, t_msmt, t_csd, size(Vm,4));

    % NOT a check, a number that must travel with the result: the two arms'
    % constraint strengths, and the l >= 2 band power that shows whether they
    % came back comparably sharp. A large gap here means the arms are not
    % being compared like for like, whatever the peak table says.
    Lv6 = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
    pw  = @(v) arrayfun(@(il) norm(v(Lv6==il))/v(1), (2:2:Lf));
    fprintf('        band power l>=2, mean over voxels:  SSST %s   MSMT %s\n', ...
            mat2str(round(pw(mean(SH_CSD{iL},1))*1e3)/1e3), ...
            mat2str(round(pw(mean(SH_MSMT{iL},1))*1e3)/1e3));
    fprintf('        (msmt_csd -neg_lambda %g -norm_lambda %g; csd constraint strength 1)\n', ...
            MSMT_NEG_LAMBDA, MSMT_NORM_LAMBDA);

    % CHECK. The right number of coefficients, and nothing non-finite: a NaN in
    % an SH volume breaks downstream tractography silently.
    ok_nc  = (size(Vm,4) == nc) && (size(Vc,4) == nc);
    ok_fin = all(isfinite(SH_MSMT{iL}(:))) && all(isfinite(SH_CSD{iL}(:)));
    fprintf('   CHECK Lmax %d: %d coefficients, all finite            %s\n', ...
            Lf, nc, VERDICT{1+(ok_nc && ok_fin)});

    c0 = [mean(SH_MSMT{iL}(:,1)), mean(Vg(:)), mean(Vf(:))];
    tissue_share(iL,:) = c0/sum(c0);
end

fprintf('   MSMT-CSD tissue share (mean l = 0 coefficient, normalised):\n');
fprintf('        Lmax        WM        GM       CSF\n');
for iL = 1:numel(LMAX_LIST)
    fprintf('        %4d  %8.4f  %8.4f  %8.4f\n', LMAX_LIST(iL), tissue_share(iL,:));
end
fprintf(['        Every voxel here is pure %s tissue with fw = %.2f, so anything\n' ...
         '        outside the WM column is leakage. Section 6.2 of the Monte Carlo\n' ...
         '        report found MSMT''s CSF fraction is NOT a usable free-water estimate.\n'], ...
        C.PRESET, K(5));

ARMS{end+1} = struct('name', 'SSST-CSD', 'sh', {SH_CSD});
ARMS{end+1} = struct('name', 'MSMT-CSD', 'sh', {SH_MSMT});
end

NARM = numel(ARMS);
fprintf('\n   %d arm(s) to score:', NARM);
for ia = 1:NARM, fprintf('  %s', ARMS{ia}.name); end
fprintf('\n\n');

%% Step 7 -- peaks, angular error and spurious peaks, per arm, Lmax and SNR
% The only questions a tractography algorithm asks of an fODF are "how many
% fibres, and pointing where", so those are what get scored -- and *every arm
% goes through these same lines*, so a difference between arms is the
% deconvolution rather than the scoring.
%
% Peaks are found by evaluating the fODF on a dense direction set and keeping
% every direction not smaller than any neighbour within |PEAK_NBR| degrees,
% then keeping those whose *anisotropic* amplitude is at least |PEAK_REL| of the
% largest. Subtracting the isotropic part matters: without it a nearly isotropic
% fODF looks like it has many strong peaks.
%
% *The isotropic part is each voxel's OWN l = 0 term, not the constant
% 1/(4*pi).* That constant is the isotropic part only for an fODF in SMI's
% |p_00 = 1| convention. An MRtrix FOD is unnormalised -- its l = 0 coefficient
% varies per voxel and carries apparent fibre density -- so subtracting the
% constant would silently mis-threshold both CSD arms. For the SMI arm the two
% are the same number, and that is checked below rather than assumed.
%
% The finder itself lives in |helpers/fODF_peak_score.m| rather than inline, so
% that this file and anything else scoring an fODF provably run the same code.
%
% Each Lmax is scored against *its own ceiling*: the ground truth truncated to
% that same Lmax, with no noise and no fitting. That is what separates "the
% method failed" from "this angular order cannot represent the answer". The
% ceiling is not scored by a copy of the peak finder -- the truth is prepended
% as the *first row of every block*, so it goes through the same code as the
% realisations do.
%
% *Four numbers per cell*, because a sweep is where they stop agreeing:
%
% * *correct count* -- the fraction of realisations recovering exactly the true
%   number of fibres. What tractography actually consumes.
% * *mean error* -- the bias: how far the primary peak sits from the nearest
%   true axis, averaged. This is the mean of a non-negative quantity, so it does
%   not go to zero even for an unbiased estimator; read it against the noise-free
%   row of the same curve, not against zero.
% * *std error* -- the spread of that same quantity across realisations.
% * *spurious* -- peaks found beyond the true count, per voxel.
%
% *One weakness worth knowing before reading the numbers.* The angular error
% uses the LARGEST peak only, and a symmetric crossing has two lobes of
% near-equal amplitude, so which one is "primary" can be decided by the last few
% ulps. The reported error can then jump between the two lobes' errors without
% the reconstruction meaningfully changing. Read it alongside the peak count and
% the spurious count, never on its own.

PEAK_NBR = 12;      % degrees, neighbourhood for the local-maximum test
PEAK_REL = 0.30;    % keep peaks at >= 30% of the largest anisotropic amplitude
PK   = fODF_peak_score();
de   = H.dirs(1500);
ND   = size(de,1);
ctxP = PK.setup(de, PEAK_NBR, PEAK_REL);

fprintf('Step 7: peaks (%d directions, %d deg neighbourhood, %.0f%% threshold)\n', ...
        ND, PEAK_NBR, 100*PEAK_REL);
fprintf('   The grid spacing sets a floor on angular error of roughly %.1f deg;\n', ...
        sqrt(4*pi/ND)*180/pi/2);
fprintf('   the ceiling line under each condition shows that floor directly.\n');

res_all  = zeros(NARM, numel(LMAX_LIST), NSNR, NCOND);  % % with the true count
bias_all = nan(NARM, numel(LMAX_LIST), NSNR, NCOND);    % mean angular error, deg
sd_all   = nan(NARM, numel(LMAX_LIST), NSNR, NCOND);    % its standard deviation
med_all  = nan(NARM, numel(LMAX_LIST), NSNR, NCOND);    % median angular error
spur_all = zeros(NARM, numel(LMAX_LIST), NSNR, NCOND);  % spurious peaks / voxel
amp_all  = nan(NARM, numel(LMAX_LIST), NSNR, NCOND);    % anisotropic peak amplitude
iso_all  = nan(NARM, numel(LMAX_LIST), NSNR, NCOND);    % the isotropic term
ceil_n   = zeros(NARM, numel(LMAX_LIST), NCOND);        % peaks the truth gives
ceil_err = nan(NARM, numel(LMAX_LIST), NCOND);          % and the truth's error

for ia = 1:NARM
    fprintf('\n   ===== %s =====\n', ARMS{ia}.name);
    for iL = 1:numel(LMAX_LIST)
        Lf = LMAX_LIST(iL);
        Ye = SMI.get_even_SH(de, Lf, C.CS_PHASE);
        nc = size(Ye,2);
        fprintf('   Lmax %d\n', Lf);
        fprintf('     condition    SNR   correct count   mean err   std err   median err   spurious\n');

        for ic = 1:NCOND
            ntrue = numel(axes_gt{ic});
            for k = 1:NSNR
                is   = SNR_ORD(k);
                rows = find(cond_id == ic & snr_id == is);

                % Row 1 is the ground truth truncated to THIS Lmax; rows 2:end
                % are the realisations of this condition at this SNR.
                blk = [sh_gt(ic,1:nc); ARMS{ia}.sh{iL}(rows,:)];
                [nfound, aerr, apk, a0] = PK.score(blk, Ye, axes_gt{ic}, ctxP);

                ceil_n(ia,iL,ic)   = nfound(1);       % identical at every SNR
                ceil_err(ia,iL,ic) = aerr(1);
                nf  = nfound(2:end);
                ae  = aerr(2:end);
                ap  = apk(2:end);
                ai  = a0(2:end);
                fin = isfinite(ae);

                res_all(ia,iL,is,ic)  = 100*mean(nf == ntrue);
                spur_all(ia,iL,is,ic) = mean(max(nf - ntrue, 0));
                if any(fin)
                    bias_all(ia,iL,is,ic) = mean(ae(fin));
                    med_all(ia,iL,is,ic)  = median(ae(fin));
                end
                if sum(fin) > 1
                    sd_all(ia,iL,is,ic) = std(ae(fin));
                end
                % Amplitude, in this arm's OWN convention. Never compare it
                % across arms -- an SMI fODF integrates to 1 and an MRtrix FOD
                % does not -- only within an arm, across tissue or SNR.
                if any(isfinite(ap)), amp_all(ia,iL,is,ic) = mean(ap(isfinite(ap))); end
                if any(isfinite(ai)), iso_all(ia,iL,is,ic) = mean(ai(isfinite(ai))); end

                fprintf('     %9s  %5s   %11.1f%%   %8.2f  %8.2f     %8.2f   %8.3f\n', ...
                        COLHDR{ic}, SNR_LABEL{is}, res_all(ia,iL,is,ic), ...
                        bias_all(ia,iL,is,ic), sd_all(ia,iL,is,ic), ...
                        med_all(ia,iL,is,ic), spur_all(ia,iL,is,ic));
            end

            if ceil_n(ia,iL,ic) == ntrue
                ctxt = 'resolvable';
            else
                ctxt = sprintf('%d of %d -- NOT resolvable', ceil_n(ia,iL,ic), ntrue);
            end
            % %12s, not %11s: the data rows above carry a literal %% after their
            % %11.1f, so the ceiling row needs one more column to stay aligned.
            fprintf('     %9s  truth   %12s   %8.2f                          %8.3f   <- ceiling, %s\n', ...
                    COLHDR{ic}, '--', ceil_err(ia,iL,ic), ...
                    max(ceil_n(ia,iL,ic)-ntrue, 0), ctxt);
        end
    end
end
fprintf('\n');
fprintf('   Reading this table: a low "correct count" next to a ceiling that says NOT\n');
fprintf('   resolvable is the angular order failing, not the method. A low count next to\n');
fprintf('   a resolvable ceiling is the method or the noise, and the SNR column says\n');
fprintf('   which -- if it is still low at SNR inf, noise was never the problem.\n');
fprintf('   The ceiling is a property of the TRUTH, so it is identical across arms;\n');
fprintf('   it is printed under each because that is where it is read.\n\n');

% ---- amplitude, which the angular table cannot show
% Tractography terminates on fODF AMPLITUDE, not on angular error, so an arm can
% score a perfect orientation and still be untrackable. This block is the only
% place that is visible.
%
% Read DOWN a column, never across: the numbers are in each arm's own
% convention. SMI's fODF has p_00 = 1 and integrates to 1; an MRtrix FOD is
% unnormalised and its amplitude carries apparent fibre density. The ratio
% between two conditions within one arm is the scale-free quantity.
fprintf('   peak amplitude above the isotropic floor, mean over voxels\n');
fprintf('   (each arm in ITS OWN scale -- compare down a column, not across)\n');
fprintf('     %-10s %-5s', 'arm', 'Lmax');
for k = 1:NSNR, fprintf('%11s', SNR_LABEL{SNR_ORD(k)}); end
fprintf('   iso term\n');
for ia = 1:NARM
    for iL = 1:numel(LMAX_LIST)
        fprintf('     %-10s %-5d', ARMS{ia}.name, LMAX_LIST(iL));
        for k = 1:NSNR
            fprintf('%11.4f', amp_all(ia,iL,SNR_ORD(k),1));
        end
        fprintf('%11.4f\n', iso_all(ia,iL,SNR_ORD(end),1));
    end
end
fprintf('\n');

% CHECK. The isotropic subtraction. For an fODF in SMI's convention the
% per-voxel l = 0 term IS the constant 1/(4*pi); for an MRtrix FOD it is not.
% Verify the general form reproduces the constant on the SMI arm, so changing
% the form to accommodate MRtrix cannot have broken the arm that was already
% right. Not required to be exactly zero: sh(:,1) is 1/sqrt(4*pi), so the
% general form evaluates (1/sqrt(4*pi))*(1/sqrt(4*pi)) where the constant form
% evaluates 1/(4*pi) -- one rounding step apart.
e_iso = max(abs(ARMS{1}.sh{1}(:,1)/sqrt(4*pi) - 1/(4*pi)));
fprintf('   CHECK SMI''s own l=0 term equals 1/(4*pi)   max|err| = %.2e   %s\n', ...
        e_iso, VERDICT{1+(e_iso < 1e-15)});

% CHECK, and a measurement: an MRtrix FOD really is on a different scale, which
% is *why* the constant would have been wrong. Measured rather than asserted.
if NARM > 1
    l0 = ARMS{2}.sh{1}(:,1);
    fprintf('   %s l=0 coefficient: min %.4f, max %.4f, std %.4f\n', ...
            ARMS{2}.name, min(l0), max(l0), std(l0));
    fprintf('   CHECK an MRtrix FOD is unnormalised (l=0 varies per voxel)  %s\n', ...
            VERDICT{1+(std(l0) > 0)});
    fprintf('        SMI''s is the fixed %.4f, which is the whole reason the peak\n', ...
            1/sqrt(4*pi));
    fprintf('        finder subtracts each voxel''s own l=0 rather than a constant.\n');
end

% CHECK. The noise-free arm must recover the true fibre count wherever the
% truth itself resolves. With no noise the only things left between the fit and
% the truth are the band limit and the estimated kernel, so a failure here makes
% nothing at a finite SNR interpretable. Conditions whose ceiling is already
% below the true count are *skipped*, not failed: there the truth does not
% resolve either, and demanding that the fit do better than the truth would be
% the wrong test.
i_inf = find(isinf(SNR_LIST), 1);
if ~isempty(i_inf)
    ok_nf = true; n_tested = 0;
    for ia = 1:NARM
        for iL = 1:numel(LMAX_LIST)
            for ic = 1:NCOND
                if ceil_n(ia,iL,ic) == numel(axes_gt{ic})
                    ok_nf = ok_nf && (res_all(ia,iL,i_inf,ic) > 99);
                    n_tested = n_tested + 1;
                end
            end
        end
    end
    fprintf('   CHECK noise-free arms recover the true count (%d of %d cells resolve)  %s\n', ...
            n_tested, NARM*numel(LMAX_LIST)*NCOND, VERDICT{1+ok_nf});

    % CHECK. Determinism of the noise-free block, tested on the SH coefficients
    % directly.
    %
    % *This check used to read |std(angular error) == 0| and that was wrong in
    % two independent ways*, both found by actually running it:
    %
    % # |std(x) == 0| is not a determinism test. |std| computes |sum(x)/n|
    %   first, and for n identical values that division need not return the
    %   value exactly, so a perfectly deterministic block can yield ~1e-16.
    %   Whether it does depends on the bit pattern of the particular angular
    %   error, which is why it could fail for one arm and pass for another on
    %   identical logic.
    % # |SMI.fit| is genuinely NOT bit-reproducible voxel to voxel. Measured:
    %   27 voxels on a [3 3 3] grid fed a bit-identical signal at sigma = 0 come
    %   back as three distinct answers, one per slice, differing by 1.6e-12 in
    %   plm. That is BLAS blocking -- reduction order depends on where a voxel
    %   sits in the array -- and it is twelve orders of magnitude below the peak
    %   finder's 2.6 degree grid resolution.
    %
    % Both MRtrix arms ARE bit identical, because dwi2fod works one voxel at a
    % time. So the check is kept rather than dropped: it separates a method that
    % is exactly reproducible voxel to voxel from one that is reproducible only
    % to floating point, and it prints the measured number either way.
    for ia = 1:NARM
        worst = 0;
        for iL = 1:numel(LMAX_LIST)
            for ic = 1:NCOND
                rows = find(cond_id == ic & snr_id == i_inf);
                blk  = ARMS{ia}.sh{iL}(rows,:);
                worst = max(worst, max(max(abs(blk - repmat(blk(1,:), numel(rows), 1)))));
            end
        end
        fprintf('   CHECK %-9s noise-free block reproducible  max|row-row1| = %.1e   %s\n', ...
                ARMS{ia}.name, worst, VERDICT{1+(worst < 1e-9)});
    end
end

% CHECK. More noise must not help -- per arm, and the exceptions are the
% interesting part. Compared at the extremes rather than pairwise so that Monte
% Carlo error between adjacent SNRs cannot trip it.
%
% *This is a check on the simulation, not a law about methods.* Its premise is
% that an arm's error is dominated by noise, so removing the noise must not make
% it worse; a swapped SNR label or a mis-indexed block would show up here. But an
% arm whose error is dominated by SYSTEMATIC bias -- a response mismatch, say --
% can genuinely do better with noise than without, because noise jitters its
% primary lobe around and sometimes lands it nearer a true axis than the biased
% answer sits. MSMT-CSD at 60 degrees with a blunt response does exactly that.
%
% So a violation is classified rather than simply failed. An arm whose noise-free
% error is already far above the ceiling is bias limited and is reported as a
% FINDING; an arm sitting near its ceiling that still gets worse without noise
% has a real problem, and that is the one that fails.
is_lo = SNR_ORD(1); is_hi = SNR_ORD(end);
ok_mono = true;
for ia = 1:NARM
    for iL = 1:numel(LMAX_LIST)
        for ic = 1:NCOND
            a_hi = bias_all(ia,iL,is_hi,ic);
            a_lo = bias_all(ia,iL,is_lo,ic);
            if ~isfinite(a_hi) || ~isfinite(a_lo), continue, end
            if a_hi <= a_lo, continue, end
            % Violated. Bias limited, or a genuine failure?
            cE = ceil_err(ia,iL,ic);
            if isfinite(cE) && a_hi > 2*cE
                fprintf(['   NOTE  %-9s Lmax %d %s: error at SNR %s is %.2f deg,\n' ...
                         '         ABOVE its %.2f deg at SNR %s and %.1fx its %.2f deg ceiling.\n' ...
                         '         Bias limited, not a simulation fault: with no noise left, what\n' ...
                         '         remains is the response mismatch, and noise partly masks it.\n'], ...
                        ARMS{ia}.name, LMAX_LIST(iL), COLHDR{ic}, ...
                        SNR_LABEL{is_hi}, a_hi, a_lo, SNR_LABEL{is_lo}, a_hi/cE, cE);
            else
                fprintf(['   %-9s Lmax %d %s: error at SNR %s (%.2f) exceeds SNR %s (%.2f)\n' ...
                         '         while sitting near its %.2f deg ceiling -- a real failure.\n'], ...
                        ARMS{ia}.name, LMAX_LIST(iL), COLHDR{ic}, ...
                        SNR_LABEL{is_hi}, a_hi, SNR_LABEL{is_lo}, a_lo, cE);
                ok_mono = false;
            end
        end
    end
end
fprintf('   CHECK no arm sitting at its ceiling gets worse without noise   %s\n\n', ...
        VERDICT{1+ok_mono});

%% Step 8 -- export the fODFs so MRtrix can check them
% Everything above is scored by this file. That is fine for reading, but it is
% not independent: the same code wrote the fODF and found its peaks. The point
% of this step is to hand the fODFs to something that has never seen this
% package.
%
% The fit ran at |CS_phase = 0|, where SMI's spherical harmonic basis *is*
% MRtrix's, so the coefficients need no conversion -- they are written straight
% out as MRtrix SH images, one per arm per Lmax.
%
% *Every arm is exported, not just SMI*, so |sh2peaks| can be run over all three
% and its answers compared against the table Step 7 printed. That comparison is
% the independent check: a basis error, which is the failure mode this is really
% guarding against, would show as ~71 degrees on every arm at once.

% The exported image holds the WHOLE sweep -- every SNR block, in the order
% Step 5 laid them out -- so a single sh2peaks run covers all noise levels and
% the key below says which voxel was which.
edir = fullfile(pkgdir, 'export');
if ~exist(edir, 'dir'), mkdir(edir); end
GRID_ALL = C.pick_grid(NVOX);
MR = mrtrix_io();
fprintf('Step 8: writing MRtrix SH images to %s\n', edir);
for ia = 1:NARM
    nm = strrep(lower(ARMS{ia}.name), '-', '');
    for iL = 1:numel(LMAX_LIST)
        Lf = LMAX_LIST(iL);
        fn = fullfile(edir, sprintf('%sfod_%s_lmax%d', nm, C.PRESET, Lf));
        MR.write([fn '.mif'], reshape(ARMS{ia}.sh{iL}, [GRID_ALL size(ARMS{ia}.sh{iL},2)]));
        fprintf('   %sfod_%s_lmax%d.mif  [%s x %d]\n', nm, C.PRESET, Lf, ...
                mat2str(GRID_ALL), size(ARMS{ia}.sh{iL},2));
    end
end

% A voxel-order key, so a peak found in MRtrix can be matched to the condition
% AND the noise level that produced it. Column 1 is the linear voxel index in
% MRtrix's order, column 2 the SNR, column 3 the crossing angle, columns 4-6 and
% 7-9 the true fibre axes.
key = zeros(NVOX, 9);
for v = 1:NVOX
    ic = cond_id(v);
    key(v,1) = v;
    key(v,2) = SNR_LIST(snr_id(v));
    key(v,3) = C.ANGLES(ic);
    key(v,4:6) = axes_gt{ic}{1};
    if numel(axes_gt{ic}) == 2, key(v,7:9) = axes_gt{ic}{2}; end
end
fid = fopen(fullfile(edir, sprintf('voxel_key_%s.txt', C.PRESET)), 'w');
fprintf(fid, '%% voxel  SNR  crossing_deg  axis1_x axis1_y axis1_z  axis2_x axis2_y axis2_z\n');
fprintf(fid, '%% axis2 is 0 0 0 for the single fibre condition. SNR Inf is the noise-free arm.\n');
fprintf(fid, '%% Voxel order is column-major over the %s grid, matching the .mif images:\n', ...
        mat2str(GRID_ALL));
fprintf(fid, '%% %d contiguous blocks of %d voxels, one per SNR, in the order %s.\n', ...
        NSNR, NVOX_SNR, strjoin(SNR_LABEL, ', '));
fprintf(fid, '%d %g %d %.9f %.9f %.9f %.9f %.9f %.9f\n', key');
fclose(fid);
fprintf('   voxel_key_%s.txt  SNR and true fibre axes per voxel, for matching peaks back\n\n', C.PRESET);

fprintf('   To check these against MRtrix3, from %s:\n', edir);
fprintf('     sh2peaks smifod_%s_lmax6.mif peaks_lmax6.mif -num 4\n', C.PRESET);
fprintf('     mrinfo   smifod_%s_lmax6.mif\n', C.PRESET);
fprintf('     mrconvert smifod_%s_lmax6.mif out.mif   # if you prefer a single file\n', C.PRESET);
fprintf('     mrview   smifod_%s_lmax6.mif -odf.load_sh smifod_%s_lmax6.mif\n', C.PRESET, C.PRESET);
fprintf('   sh2peaks writes 3 components per peak; compare against voxel_key.txt.\n');
fprintf('   Nothing in this walkthrough runs those commands -- that is the point.\n\n');

%% End of the per-kernel loop
% Stash everything the figures need. Kernel-independent quantities (the ground
% truth, the protocol, the shell assignment) are computed once above and are not
% duplicated here.

RUN{ik} = struct( ...
    'name',    C.PRESET, ...
    'K',       C.K_WM, ...
    'S_clean', S_clean, ...
    'S_noisy', S_noisy, ...
    'fits',    {fits}, ...
    'arms',    {ARMS}, ...            % SMI first, then the CSD arms
    'armname', {cellfun(@(a) a.name, ARMS, 'UniformOutput', false)}, ...
    'res',     res_all, ...           % all now [NARM x NLMAX x NSNR x NCOND]
    'bias',    bias_all, ...
    'sd',      sd_all, ...
    'spur',    spur_all, ...
    'amp',     amp_all, ...
    'iso',     iso_all, ...
    'ceil_n',  ceil_n);

end   % ik, the kernel loop
fprintf('==================== both kernels done ====================\n\n');

%% Figures
% Six figures, plus one more per kernel comparing the arms. Set
% |MAKE_FIGURES = false| to skip them.
%
% All glyphs use |SMI_response_helpers|, the renderer MRtrix's |shview| logic
% implies: *radius is |amplitude|, colour is the signed amplitude*, so negative
% lobes show as a colour change instead of being folded silently into the
% surface. The band-limited truth really is negative over much of the sphere
% (Step 2), and a renderer that hid that would hide the most surprising result
% in the file.
%
% Every 3D panel opens *isometric*, so all subfigures are readable without
% touching the camera.
%
%  Fig 1  the ground truth fODF, then the response glyph of each kernel per
%         shell -- what is being convolved, next to what it is convolved with
%  Fig 2  the spherical signal, HEALTHY: b down the rows, SNR across
%  Fig 3  the spherical signal, EDEMA: same layout, same radial scale
%  Fig 4  reconstructed fODFs, HEALTHY: Lmax down, truth then one column per SNR
%  Fig 5  reconstructed fODFs, EDEMA: same layout
%  Fig 6  bias, spread and spurious peaks against SNR -- one ROW PER KERNEL

MAKE_FIGURES = true;
ISO_VIEW     = [1 1 1];      % isometric camera: equal foreshortening on x, y, z
GLYPH_N      = 121;          % mesh resolution for every glyph
GLYPH_NEG    = 'clamp';      % what a glyph does with NEGATIVE fODF amplitude.
                             %   'clamp'  radius = max(amplitude,0)
                             %   'abs'    radius = |amplitude|, the shview convention
                             %
                             % A band-limited fODF rings and the rings go below
                             % zero. Measured on a noise-free edema voxel at
                             % Lmax 6, MSMT-CSD's fODF is negative over 64.6% of
                             % the sphere and SSST-CSD's over 48.6%. Under 'abs'
                             % every one of those regions is drawn at POSITIVE
                             % radius, so the glyph grows lobes that are not
                             % fibres and reads as inflated -- and the sign is
                             % only in the colour, which flat-shaded gnuplot
                             % makes easy to miss. 'clamp' shows the fODF a
                             % tractography algorithm would actually follow.
LMAX_FIG     = 6;            % the angular order Figure 7 compares the arms at.
                             % 6 because that is what the published comparison
                             % in Reports/ ran every arm at, so Figure 7 is the
                             % row to read against the report. Falls back to the
                             % first entry of LMAX_LIST if 6 is not in it.
GLYPH_LIM    = [-1.02 1.02]; % axis limits, applied to FIGURE 1 ONLY

% *Why GLYPH_LIM exists, and why it is not optional where it is used.*
% |axis equal| fixes the ASPECT RATIO of a panel; it does not fix the axis
% LIMITS. Left to itself MATLAB autoscales each subplot to its own data, so a
% glyph half the size of its neighbour gets an axis range half as wide and ends
% up drawn exactly the same size on the page. Scaling the radius then has no
% visible effect at all. Both halves are needed: scale the radius AND pin the
% limits.
%
% *This is deliberately applied to Figure 1 alone.* Figure 1 is the one panel
% whose point is that the edema response is genuinely smaller than the healthy
% one, so there size must carry meaning. Everywhere else a shared scale would
% just shrink one figure into illegibility for no gain: Figures 2 and 3 are
% about the SHAPE of the signal at each shell and SNR, and Figures 4 and 5
% about the shape of the reconstruction, so every panel there is left to
% autoscale and fill its box as it did before.

if MAKE_FIGURES
    fprintf('Figures: rendering ...\n');
    [THg, PHg, dirs_g] = RH.grid(GLYPH_N, GLYPH_N);
    Yg   = SMI.get_even_SH(dirs_g, LMAX_GT, C.CS_PHASE);
    nsh  = numel(b_shell);
    ic   = 1;                                  % the only condition: 60 degrees

    % ================================================================
    % FIGURE 1 -- the ground truth and each kernel's response per shell
    % ================================================================
    % The left column is the fODF being convolved. It is the same for both
    % kernels: the ground truth geometry does not depend on the tissue.
    %
    % To its right, one row per kernel and one glyph per shell.
    %
    % *There is deliberately no difference row.* An earlier version carried one,
    % showing R_healthy - R_edema per shell. It is gone because the shared
    % radial scale already makes that comparison directly -- the edema glyphs
    % are visibly smaller than the healthy ones in the same units, which is the
    % claim -- while a residual glyph restated it in a form that is harder to
    % read and invited the misreading that the residual is itself a response.
    %
    % *Every glyph on this figure shares ONE radial scale*, so size carries
    % meaning: a smaller glyph is a genuinely smaller signal. That is the whole
    % point of the panel. The edema kernel attenuates faster and is far less
    % anisotropic, and both show up as shrinkage rather than only as a change
    % of shape.
    %
    % Two things follow from the shared scale and are worth expecting:
    %
    % * *The b = 0 glyphs are identical unit spheres in both rows.*
    %   r_0(0) = sqrt(4*pi) for ANY kernel, because S(0)/S0 = 1 exactly. That identity holds here because
    %   |B0_SNAP| set the acquired b = 5 s/mm^2 volumes to exactly 0; with
    %   |B0_SNAP = 0| they would differ slightly instead. Either way it is
    %   correct, not a bug.
    % * The scale is set by that b = 0 sphere, the largest thing on the figure.
    %
    % This is also the panel built to take CSD and MSMT-CSD: dwi2response
    % writes zonal coefficients in exactly the form RH.zh_glyph draws, so an
    % estimated response becomes another row here with no conversion.
    th_prof = linspace(0, pi, 361);
    nrow1   = NKERN;

    % one radial scale for the whole figure, over every glyph it will draw
    r_sh = cell(NKERN, nsh);
    rmax = 0;
    for ikk = 1:NKERN
        for i = 1:nsh
            r_sh{ikk,i} = RH.zh(KERNELS{ikk}.K, b_shell(i), LMAX_GT, C.D_FW);
            rmax = max(rmax, max(abs(RH.profile(r_sh{ikk,i}, th_prof))));
        end
    end
    if rmax <= 0, rmax = 1; end

    figure('Name', 'Fig 1  ground truth fODF and the kernel responses per shell');
    col1 = 1 + (0:nrow1-1)*(1+nsh);        % the truth spans the first column
    subplot(nrow1, 1+nsh, col1);
    % The ground truth is an fODF, not a signal, so it is NOT on the response
    % scale -- the two are different quantities and a shared scale between them
    % would mean nothing. It is normalised to its own peak instead.
    gt_pk = max(abs(sh_gt(ic,:) * Yg'));
    if gt_pk <= 0, gt_pk = 1; end
    [X, Y, Z, Cc] = RH.sh_glyph(plm_gt(ic,:), LMAX_GT, C.CS_PHASE, ...
                                GLYPH_N, GLYPH_N, 1/gt_pk, GLYPH_NEG);
    surf(X, Y, Z, Cc); shading interp;
    axis equal off vis3d; view(ISO_VIEW);
    set(gca, 'XLim', GLYPH_LIM, 'YLim', GLYPH_LIM, 'ZLim', GLYPH_LIM);
    camlight headlight; lighting gouraud;
    title(sprintf('ground truth\n%g deg, kappa %g', ANGLES(ic), C.KAPPA));

    for ikk = 1:NKERN
        for i = 1:nsh
            subplot(nrow1, 1+nsh, (ikk-1)*(1+nsh) + 1 + i);
            [X, Y, Z, Cc] = RH.zh_glyph(r_sh{ikk,i}, GLYPH_N, GLYPH_N, 1/rmax, GLYPH_NEG);
            surf(X, Y, Z, Cc); shading interp;
            axis equal off vis3d; view(ISO_VIEW);
            set(gca, 'XLim', GLYPH_LIM, 'YLim', GLYPH_LIM, 'ZLim', GLYPH_LIM);
            camlight headlight; lighting gouraud;
            title(sprintf('%s, b = %.0f', KERNELS{ikk}.name, b_shell(i)));
        end
    end

    % The same comparison as numbers, since a glyph is hard to read off a page.
    fprintf('   Fig 1: response peak amplitude on the shared scale (max = %.3f)\n', rmax);
    fprintf('        b   ');
    for ikk = 1:NKERN, fprintf('%12s', KERNELS{ikk}.name); end
    fprintf('\n');
    for i = 1:nsh
        fprintf('     %4.2f   ', b_shell(i));
        for ikk = 1:NKERN
            fprintf('%12.4f', max(abs(RH.profile(r_sh{ikk,i}, th_prof))));
        end
        fprintf('\n');
    end
    fprintf(['   Every kernel gives the same response at b = 0, since S(0)/S0 = 1, so\n' ...
             '   that row of the difference column is zero to machine precision. It is\n' ...
             '   non-zero only if B0_SNAP is disabled and the acquired b = 5 s/mm^2 is\n' ...
             '   used as-is.\n\n']);

    % ================================================================
    % FIGURES 2 and 3 -- the spherical signal, one figure per kernel
    % ================================================================
    % Radius is S(u)/S0. The noise-free column is the model evaluated on a dense
    % grid; every other column is the spherical harmonic fit of ONE
    % representative realisation at that SNR, which is the same projection SMI
    % does internally, so what is drawn is what the fit actually sees.
    %
    % Both figures share ONE radial scale, computed across both kernels, so the
    % edema signal can be compared against the healthy one by eye rather than
    % only by number.
    sig = cell(NKERN, nsh-1, NSNR);
    smax = 0;
    for ikk = 1:NKERN
        for j = 2:nsh
            for kk = 1:NSNR
                is = SNR_ORD(kk);
                if isinf(SNR_LIST(is))
                    bq = b_shell(j)*ones(1, size(dirs_g,1));
                    sg = H.signal(plm_gt(ic,:), [KERNELS{ikk}.K 1 1], bq, ...
                                  ones(size(bq)), zeros(size(bq)), dirs_g, ...
                                  LMAX_GT, C.CS_PHASE, C.D_FW);
                    A = reshape(sg, size(THg));
                else
                    v  = find(snr_id(:) == is, 1);
                    m  = (shell_id(:)' == j);
                    Ym = SMI.get_even_SH(bvecs(m,:), LMAX_GT, C.CS_PHASE);
                    cm = Ym \ RUN{ikk}.S_noisy(v, m)';
                    A  = reshape(Yg*cm, size(THg));
                end
                sig{ikk, j-1, kk} = A;
                smax = max(smax, max(abs(A(:))));
            end
        end
    end
    for ikk = 1:NKERN
        figure('Name', sprintf('Fig %d  spherical signal, %s', ...
                               1+ikk, RUN{ikk}.name));
        for j = 2:nsh
            for kk = 1:NSNR
                subplot(nsh-1, NSNR, (j-2)*NSNR + kk);
                [X, Y, Z, Cc] = RH.glyph(sig{ikk, j-1, kk}, THg, PHg, 1/smax, GLYPH_NEG);
                surf(X, Y, Z, Cc); shading interp;
                axis equal off vis3d; view(ISO_VIEW);
                camlight headlight; lighting gouraud;
                title(sprintf('b %.0f, SNR %s', b_shell(j), ...
                              SNR_LABEL{SNR_ORD(kk)}));
            end
        end
    end

    % ================================================================
    % FIGURES 4 and 5 -- the SMI reconstruction, one figure per kernel
    % ================================================================
    % *The SMI arm only*, Lmax down the rows: these two figures are about how
    % angular order and noise affect SMI's own reconstruction. The first column
    % is the band-limited ground truth at that Lmax -- the bound nothing below
    % it can beat -- then one column per SNR, worst first. Same renderer and
    % camera as Figure 1, so truth and reconstruction can be compared directly.
    % Figure 7 is where the ARMS are put side by side.
    for ikk = 1:NKERN
        figure('Name', sprintf('Fig %d  reconstructed fODFs, %s', ...
                               3+ikk, RUN{ikk}.name));
        for iL = 1:numel(LMAX_LIST)
            Lf = LMAX_LIST(iL);
            nc = (Lf/2+1)*(Lf+1);
            Lv = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
            sc = sqrt((2*Lv(2:end)'+1)/(4*pi));

            subplot(numel(LMAX_LIST), NSNR+1, (iL-1)*(NSNR+1) + 1);
            [X, Y, Z, Cc] = RH.sh_glyph(sh_gt(ic,2:nc)./sc, Lf, C.CS_PHASE, ...
                                        GLYPH_N, GLYPH_N, 1, GLYPH_NEG);
            surf(X, Y, Z, Cc); shading interp;
            axis equal off vis3d; view(ISO_VIEW);
            camlight headlight; lighting gouraud;
            title(sprintf('Lmax %d\ntruth', Lf));

            for kk = 1:NSNR
                is = SNR_ORD(kk);
                subplot(numel(LMAX_LIST), NSNR+1, (iL-1)*(NSNR+1) + 1 + kk);
                v  = find(snr_id(:) == is, 1);
                [X, Y, Z, Cc] = RH.sh_glyph(RUN{ikk}.fits{iL}.sh(v, 2:nc)./sc, ...
                                            Lf, C.CS_PHASE, GLYPH_N, GLYPH_N, 1, GLYPH_NEG);
                surf(X, Y, Z, Cc); shading interp;
                axis equal off vis3d; view(ISO_VIEW);
                camlight headlight; lighting gouraud;
                title(sprintf('SNR %s', SNR_LABEL{is}));
            end
        end
    end

    % ================================================================
    % FIGURE 6 -- the comparison figure: arms against SNR
    % ================================================================
    % *ONE FIGURE PER KERNEL. Rows are Lmax, columns are the metric, and the
    % three lines in every panel are the arms.* That is the layout the
    % comparison wants: within a panel the only thing changing is the method,
    % and the arms sit on the same axes because they were run on the same
    % voxels -- there is no Monte Carlo error between them to explain away.
    %
    % *Y limits are shared per column across BOTH kernels*, so a panel from the
    % edema figure can be read against the same panel of the healthy one
    % without rescaling by eye. That is the whole point of drawing them
    % separately rather than autoscaling each.
    %
    % The fourth column is AMPLITUDE, and it is the one that does not obey the
    % rule above: amplitude is in each arm's own convention -- SMI's fODF
    % integrates to 1, an MRtrix FOD does not -- so the three lines in that
    % panel are NOT on a common scale and their vertical ordering is
    % meaningless. What is meaningful is each line's own shape across SNR, and
    % the same line between the healthy and edema figures. The panel is titled
    % to say so.
    %
    % Drawn from exactly the arrays Step 7 printed, so tables and plots cannot
    % disagree.
    xs   = 1:NSNR;                       % SNR_ORD order: worst first, Inf last
    mets = {'bias', 'sd', 'spur', 'amp'};
    % Column titles are kept SHORT deliberately. The descriptive versions --
    % "mean angular error", "std of angular error", "spurious peaks per voxel",
    % "peak amplitude above the isotropic floor" -- overrun the panel width and
    % collide with their neighbours under gnuplot, which is the same failure the
    % figure titles hit once before. The long names live in notebooks/README.md.
    labs = {'bias (deg)', 'spread (deg)', 'spurious/vox', 'amp (own scale)'};
    NA   = numel(RUN{1}.armname);
    NM   = numel(mets);

    % Shared y limits per metric, computed over every kernel, arm and Lmax.
    ylim_all = cell(1, NM);
    for im = 1:NM
        hi = 0;
        for ikk = 1:NKERN
            M = RUN{ikk}.(mets{im});
            v = M(:,:,SNR_ORD,ic);
            v = v(isfinite(v));
            if ~isempty(v), hi = max(hi, max(v)); end
        end
        if ~isfinite(hi) || hi <= 0, hi = 1; end
        ylim_all{im} = [0 1.05*hi];
    end

    for ikk = 1:NKERN
        figure('Name', sprintf('Fig 6  arms against SNR, %s', RUN{ikk}.name));
        for iL = 1:numel(LMAX_LIST)
            for im = 1:NM
                subplot(numel(LMAX_LIST), NM, (iL-1)*NM + im);
                hold on;
                M = RUN{ikk}.(mets{im});
                for ia = 1:NA
                    plot(xs, squeeze(M(ia, iL, SNR_ORD, ic)), '-o', 'LineWidth', 1.5);
                end
                set(gca, 'XTick', xs, 'XTickLabel', SNR_LABEL(SNR_ORD));
                xlim([0.5 NSNR+0.5]); ylim(ylim_all{im}); grid on;
                if im == 1
                    ylabel(sprintf('Lmax %d', LMAX_LIST(iL)));
                end
                if iL == 1
                    title(labs{im});
                end
                if iL == numel(LMAX_LIST), xlabel('SNR'); end
                if iL == 1 && im == NM
                    legend(RUN{ikk}.armname, 'Location', 'best');
                end
            end
        end
    end

    % ================================================================
    % FIGURE 7 -- the arms side by side, one figure per kernel
    % ================================================================
    % *ONE ROW PER ARM*, at |LMAX_FIG|: the band-limited truth first, then one
    % column per SNR. This is the direct visual comparison -- same voxel, same
    % noise, three deconvolutions.
    %
    % Each MRtrix FOD is put onto SMI's |p_00 = 1| convention before the shared
    % renderer sees it, by dividing out its own l = 0 term. That is a change of
    % SCALE only: it multiplies every band by the same factor, so it preserves
    % peak orientation exactly. Without it the SMI glyphs would vanish next to
    % the CSD ones, which live on an unnormalised amplitude scale.
    iLf = find(LMAX_LIST == LMAX_FIG, 1);
    if isempty(iLf), iLf = 1; end
    Lf  = LMAX_LIST(iLf);
    nc  = (Lf/2+1)*(Lf+1);
    Lv  = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
    sc  = sqrt((2*Lv(2:end)'+1)/(4*pi));
    for ikk = 1:NKERN
        figure('Name', sprintf('Fig 7  arms side by side, %s, Lmax %d', ...
                               RUN{ikk}.name, Lf));
        for ia = 1:NA
            subplot(NA, NSNR+1, (ia-1)*(NSNR+1) + 1);
            [X, Y, Z, Cc] = RH.sh_glyph(sh_gt(ic,2:nc)./sc, Lf, C.CS_PHASE, ...
                                        GLYPH_N, GLYPH_N, 1, GLYPH_NEG);
            surf(X, Y, Z, Cc); shading interp;
            axis equal off vis3d; view(ISO_VIEW);
            camlight headlight; lighting gouraud;
            title(sprintf('%s\ntruth', RUN{ikk}.armname{ia}));

            for kk = 1:NSNR
                is = SNR_ORD(kk);
                v  = find(snr_id(:) == is & cond_id(:) == ic, 1);
                m  = RUN{ikk}.arms{ia}.sh{iLf}(v, 1:nc);
                % onto p_00 = 1, then strip the l = 0 term the renderer implies
                plm = (m(2:end) * ((1/sqrt(4*pi))/m(1))) ./ sc;
                subplot(NA, NSNR+1, (ia-1)*(NSNR+1) + 1 + kk);
                [X, Y, Z, Cc] = RH.sh_glyph(plm, Lf, C.CS_PHASE, ...
                                            GLYPH_N, GLYPH_N, 1, GLYPH_NEG);
                surf(X, Y, Z, Cc); shading interp;
                axis equal off vis3d; view(ISO_VIEW);
                camlight headlight; lighting gouraud;
                title(sprintf('SNR %s', SNR_LABEL{is}));
            end
        end
    end

    fprintf('   %d figures drawn, all 3D panels opening isometric.\n', 6 + NKERN);
    fprintf('   Radius is |amplitude|, colour is the SIGNED amplitude.\n\n');
end

%% Step 9 -- from here to a full campaign
% Everything above is the SMI arm alone, at |NREP| realisations across
% |numel(SNR_LIST)| noise levels. The published campaign is the same forward
% model at |NREP = 10000|, and it also runs CSD and MSMT-CSD, driven by
% |gen_montecarlo.m| and |run_mrtrix.sh| rather than by hand -- one invocation
% per SNR, which is the same decomposition Step 6 uses:
%
%  cd deconv_comparison
%  octave --eval "run('oct_path.m'); gen_montecarlo(50, 10000, 'snr50')"
%  ./run_mrtrix.sh fit snr50            % CSD, MSMT-CSD and every peak, in MRtrix
%  python3 tables.py nf snr50:50 ...    % the tables in Reports/deconv_tables.md
%
% *The published tables predate the real protocol.* Everything in |Reports/| was
% measured on the older synthetic scheme, which is still on disk at
% |protocol/hcp_like_3shell.txt| and marked superseded. |gen_montecarlo.m| and
% |sweep_nonneg.m| now read the same real HCP protocol this walkthrough does, so
% regenerating those tables is a re-run rather than a code change.
%
% *Where CSD and MSMT-CSD go.* Figure 2 is the panel built to take them: both
% are driven by |run_mrtrix.sh|, and |dwi2response| writes zonal coefficients in
% exactly the form |RH.zh_glyph| already draws, so each response becomes another
% row of that montage with no conversion. Their fODFs are SH images in MRtrix's
% basis, which at |CS_phase = 0| is SMI's basis, so they drop into Figures 4 and
% 5 as further columns, and into Figure 6 as further curves.
%
% *What to change first if you want to probe it.* Everything is in the
% Configuration block at the top of this file, and nothing outside it needs
% editing: |NREP| sets the runtime, |SNR_LIST| and |LMAX_LIST| the sweep,
% |KERNELS| the tissues simulated, |KAPPA| the fibre dispersion, |ANGLES| the
% crossing, |REG| the regularizer, |PROTOCOL| the acquisition.
%
% Those settings are *local to this notebook*. Changing them here does not
% change |gen_montecarlo.m| or |sweep_nonneg.m|, which keep their own values in
% |mc_config.m| -- so a manuscript figure and the published pipeline can be
% retuned independently. What the two still share is the geometry and protocol
% *utilities*, which is what stops the fibre-axis convention and the acquisition
% reader from drifting apart.

fprintf('=== walkthrough complete ===\n');
fprintf('%d realisations per condition per SNR, %d SNR levels, Lmax %s\n\n', ...
        NREP, NSNR, mat2str(LMAX_LIST));

% One block per metric, per condition, per arm. These arrays carry a leading ARM
% index since Step 6b, so every read of them needs four subscripts -- this table
% is the last place that was still indexing them as (Lmax, SNR, condition), and
% it only shows up once the tables have already printed.
sumnames = {'correct fibre count (%)', ...
            'bias: mean angular error (deg)', ...
            'std of angular error (deg)', ...
            'spurious peaks per voxel'};
sums = {res_all, bias_all, sd_all, spur_all};
armname = RUN{NKERN}.armname;
for im = 1:numel(sums)
    fprintf('%s\n', sumnames{im});
    for ic = 1:NCOND
        fprintf('  %s\n', COLHDR{ic});
        fprintf('   %-10s %-6s', 'arm', 'Lmax');
        for k = 1:NSNR, fprintf('%10s', SNR_LABEL{SNR_ORD(k)}); end
        fprintf('\n');
        for ia = 1:numel(armname)
            for iL = 1:numel(LMAX_LIST)
                fprintf('   %-10s %-6d', armname{ia}, LMAX_LIST(iL));
                for k = 1:NSNR
                    fprintf('%10.3f', sums{im}(ia, iL, SNR_ORD(k), ic));
                end
                fprintf('\n');
            end
        end
    end
    fprintf('\n');
end
fprintf(['(this last table is the FINAL kernel only, %s -- the per-kernel\n' ...
         ' tables are printed by Step 7 inside the loop, and Figures 6 and 7\n' ...
         ' carry both.)\n\n'], RUN{NKERN}.name);
