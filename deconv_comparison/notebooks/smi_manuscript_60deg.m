%% Manuscript figures: SMI on single fibres and 60 degree crossings
% The manuscript version of |smi_simulation_walkthrough.m|, cut down to the two
% configurations the manuscript uses -- *a single fibre and a 60 degree
% crossing* -- and swept across *a list of SNRs* rather than run at one.
% Everything else -- forward model, checks, fit, export -- is the same code and
% the same |mc_config.m| constants.
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
% method.
%
% *This file is meant to be read while it runs.* It is a plain |.m| script, so
% it runs as-is in MATLAB and in GNU Octave; MATLAB can also open it and save
% it as a Live Script (|.mlx|) with no edits, since the sections and markup
% below are publish-style.
%
% *Nothing here is a reimplementation.* The forward model comes from
% |helpers/fODF_modulation_helpers.m|, the fit is the real |SMI.fit|, and the
% experiment's constants come from |mc_config.m| -- the same file
% |gen_montecarlo.m| reads. If a number here disagrees with the pipeline, the
% pipeline is what is wrong.
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
%  Step 7  peaks, angular error, spurious   -> where does it stop working?
%  Step 8  export the fODFs for MRtrix      -> check them outside this file
%  Fig 1   ground truth fibre geometry      -> the two configurations
%  Fig 2   the kernel as a response         -> one glyph per shell
%  Fig 3   the signal                       -> b down the rows, SNR across
%  Fig 4   single fibre fODFs               -> truth, then SNR across
%  Fig 5   60 degree crossing fODFs         -> truth, then SNR across
%  Fig 6   bias, std and spurious vs SNR    -> the preliminary summary plots
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

%% Setup
% |oct_path.m| puts the repository root, |helpers/| and this package on the
% path, and under Octave also loads the compatibility shims. It is the only
% path manipulation in this file.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
pkgdir = fileparts(here);                    % deconv_comparison/
if exist('OCTAVE_VERSION', 'builtin')        % quieten Octave's package-load noise
    warning('off', 'all');                   % so the CHECK lines are easy to find
end
run(fullfile(pkgdir, 'oct_path.m'));

NREP      = 1000;               % realisations per condition PER SNR. Report used 10000.
SNR_LIST  = [5 10 20 30 50 Inf];% 1/sigma at each noise level, in any order
LMAX_LIST = [4 6 8];            % angular orders to fit at -- see the next section
                                % the protocol itself is named in mc_config.m

% |Inf| is the noise-free arm and is carried as a genuine infinity rather than
% as a large finite number: sigma = 1/Inf = 0, the noise draw below reduces to
% the signal exactly, and SMI's Rician bias correction subtracts zero. The
% pipeline's own noise-free arm instead uses SNR = 1e4
% (|gen_montecarlo(1e4, 500, 'nf')|), which is the same thing to within the
% band-limiting error. The only cost of a true zero is that anything dividing
% by sigma has to be special-cased, and the checks below do that explicitly
% rather than printing NaN.
%
% Two things follow from |Inf| being exactly noiseless, and they are worth
% knowing before you wait for the run: every realisation in that block is the
% same signal, so all |NREP| voxels come back identical, and the block therefore
% costs a full fit to produce one distinct answer. It earns its place as the
% reference column of Figures 3 to 5 and as the check that the fit reaches the
% truth at all -- but if runtime is what is hurting, it is the cheapest entry to
% drop, or to keep while lowering |NREP|.
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

C = mc_config();                % kernel, dispersion, seeds, protocol
H = fODF_modulation_helpers();

% The ONLY thing this file overrides in mc_config is the condition list: the
% manuscript uses a single fibre and a 60 degree crossing. The fibre axes still
% come from C.AXIS1 and C.rotate_about, so the geometry convention is shared
% with the pipeline and cannot drift.
ANG = [0 60];
AX  = cell(1, numel(ANG));
for ic = 1:numel(ANG)
    if ANG(ic) == 0
        AX{ic} = {C.AXIS1};
    else
        AX{ic} = {C.AXIS1, C.rotate_about(C.AXIS1, ANG(ic))};
    end
end
C.ANGLES = ANG;
C.condition_axes = @(ic) AX{ic};
VERDICT = {'** FAILED **', 'ok'};             % VERDICT{1+condition}

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
% The published comparison in |Reports/| ran every arm at *Lmax 6*
% (|dwi2fod ... -lmax 6| for CSD and MSMT-CSD), so 6 is the row to compare
% against the report. 4 is what the real-data driver |run_smi_batch_mod.m|
% uses. 8 is the report's own upper control.

LMAX_GT = C.LMAX_GT;
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
% * *The b = 0 volumes are not b = 0.* They are b = 5 s/mm^2, and they carry
%   unit direction vectors like every other volume even though the direction is
%   meaningless there.
% * *The b values jitter within each shell* -- 18 distinct values across the
%   scheme. The forward model below uses the exact per-volume b; SMI bins them
%   into shells for the kernel fit, and Step 1 checks that it bins them the way
%   a human would.

% Read through mc_config, which every other arm also uses, so there is one
% definition of what was acquired. It prints a loud warning if the .bvec is not
% unit -- expect one here, and see the note under Step 3 for why it matters.
[bvals, bvecs] = C.load_protocol();
Ndwi = numel(bvals);

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
% mc_config normalised them on load -- the number that matters is the one in the
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

K  = C.K_WM;
RH = SMI_response_helpers();
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
fprintf('   CHECK S(b~0) is within 1%% of 1   max|S-1| = %.2e   %s\n', ...
        e_s0, VERDICT{1+(e_s0 < 0.01)});

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
% *Each SNR is seeded separately*, from |mc_config.m|'s seed offset by the SNR
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
for is = 1:NSNR
    rows  = find(snr_id == is);
    rblk  = S_noisy(rows,:) - S_rep(rows,:);
    s_hat = std(rblk(:));
    sg    = SIGMA_LIST(is);
    if sg == 0
        ok_is = (max(abs(rblk(:))) == 0);
        fprintf('   CHECK SNR %-4s  noise free, residual is exactly zero        %s\n', ...
                SNR_LABEL{is}, VERDICT{1+ok_is});
    else
        rel   = abs(s_hat - sg)/sg;
        ok_is = (rel < 0.10);
        fprintf('   CHECK SNR %-4s  recovered sigma %.5f vs %.5f, %4.1f%% off  %s\n', ...
                SNR_LABEL{is}, s_hat, sg, 100*rel, VERDICT{1+ok_is});
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
        options.fODF_regularization = struct('flag_nonneg', 1, 'lambda_tikhonov', 0.3);

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

%% Step 7 -- peaks, angular error and spurious peaks, per Lmax and per SNR
% The only questions a tractography algorithm asks of an fODF are "how many
% fibres, and pointing where", so those are what get scored.
%
% Peaks are found by evaluating the fODF on a dense direction set and keeping
% every direction not smaller than any neighbour within |PEAK_NBR| degrees,
% then keeping those whose *anisotropic* amplitude -- amplitude minus the
% 1/(4*pi) floor -- is at least |PEAK_REL| of the largest. Subtracting the
% floor matters: without it a nearly isotropic fODF looks like it has many
% strong peaks.
%
% Each Lmax is scored against *its own ceiling*: the ground truth truncated to
% that same Lmax, with no noise and no fitting. That is what separates "the
% method failed" from "this angular order cannot represent the answer". The
% ceiling is not scored by a copy of the peak finder -- the truth is prepended
% as the *first row of every block*, so it goes through byte for byte the same
% code as the realisations do.
%
% *Four numbers per cell*, because a sweep is where they stop agreeing:
%
% * *correct count* -- the fraction of realisations recovering exactly the true
%   number of fibres. What tractography actually consumes.
% * *mean error* -- the bias: how far the primary peak sits from the nearest
%   true axis, averaged. Note this is the mean of a non-negative quantity, so it
%   does not go to zero even for an unbiased estimator; the noise-free row is
%   what it looks like with nothing but the direction grid and the band limit.
% * *std error* -- the spread of that same quantity across realisations.
% * *spurious* -- peaks found beyond the true count, per voxel.

PEAK_NBR = 12;      % degrees, neighbourhood for the local-maximum test
PEAK_REL = 0.30;    % keep peaks at >= 30% of the largest anisotropic amplitude
de   = H.dirs(1500);
ND   = size(de,1);
cosN = cosd(PEAK_NBR);

% The neighbour set of every direction, including itself, precomputed once.
%
% This replaces the |max(a .* nbr, [], 1)| form, which built a full ND x ND
% array for every voxel and at NREP = 1000 would dominate the runtime of the
% whole file. The two are *equivalent*, not merely similar: with |nbr| logical,
% every non-neighbour contributes exactly 0 to that maximum, so the old
% threshold is |max(0, neighbourhood max)|. It is only ever compared against
% amplitudes that have already passed |a > 0|, and for those |a >= max(0, m)|
% and |a >= m| select the same directions.
NBIDX = cell(1, ND);
for j = 1:ND
    NBIDX{j} = find((de*de(j,:)') > cosN);    % includes j itself
end

fprintf('Step 7: peaks (%d directions, %d deg neighbourhood, %.0f%% threshold)\n', ...
        ND, PEAK_NBR, 100*PEAK_REL);
fprintf('   The grid spacing sets a floor on angular error of roughly %.1f deg;\n', ...
        sqrt(4*pi/ND)*180/pi/2);
fprintf('   the ceiling line under each condition shows that floor directly.\n');

res_all  = zeros(numel(LMAX_LIST), NSNR, NCOND);  % % recovering the true count
bias_all = nan(numel(LMAX_LIST), NSNR, NCOND);    % mean angular error, deg
sd_all   = nan(numel(LMAX_LIST), NSNR, NCOND);    % std of angular error, deg
med_all  = nan(numel(LMAX_LIST), NSNR, NCOND);    % median angular error, deg
spur_all = zeros(numel(LMAX_LIST), NSNR, NCOND);  % mean spurious peaks per voxel
ceil_n   = zeros(numel(LMAX_LIST), NCOND);        % peaks the truth itself gives
ceil_err = nan(numel(LMAX_LIST), NCOND);          % and the truth's own error

for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);
    Ye = SMI.get_even_SH(de, Lf, C.CS_PHASE);
    nc = size(Ye,2);
    fprintf('\n   Lmax %d\n', Lf);
    fprintf('     condition    SNR   correct count   mean err   std err   median err   spurious\n');

    for ic = 1:NCOND
        ntrue = numel(axes_gt{ic});
        for k = 1:NSNR
            is   = SNR_ORD(k);
            rows = find(cond_id == ic & snr_id == is);

            % Row 1 is the ground truth truncated to THIS Lmax; rows 2:end are
            % the realisations of this condition at this SNR.
            A    = [sh_gt(ic,1:nc); fits{iL}.sh(rows,:)] * Ye' - 1/(4*pi);
            nrow = size(A,1);

            Amax = zeros(nrow, ND);
            for j = 1:ND
                Amax(:,j) = max(A(:, NBIDX{j}), [], 2);
            end
            ismax = (A > 0) & (A >= Amax);

            nfound = zeros(nrow,1);
            aerr   = nan(nrow,1);
            for r = 1:nrow
                lm = find(ismax(r,:)); lm = lm(:);
                if isempty(lm), continue; end
                a  = A(r,:)';
                lm = lm(a(lm) >= PEAK_REL*max(a(lm)));
                [~, o] = sort(a(lm), 'descend'); lm = lm(o);
                P = de(lm,:);
                sel = true(size(P,1),1);
                for i = 1:size(P,1)
                    if ~sel(i), continue; end
                    dup = abs(P*P(i,:)') > cosN; dup(i) = false; sel(dup) = false;
                end
                P = P(sel,:);
                nfound(r) = size(P,1);
                d = zeros(1,ntrue);
                for kk = 1:ntrue
                    d(kk) = acosd(min(abs(P(1,:)*axes_gt{ic}{kk}'), 1));
                end
                aerr(r) = min(d);
            end

            ceil_n(iL,ic)   = nfound(1);          % identical at every SNR
            ceil_err(iL,ic) = aerr(1);
            nf = nfound(2:end);
            ae = aerr(2:end);
            fin = isfinite(ae);

            res_all(iL,is,ic)  = 100*mean(nf == ntrue);
            spur_all(iL,is,ic) = mean(max(nf - ntrue, 0));
            if any(fin)
                bias_all(iL,is,ic) = mean(ae(fin));
                med_all(iL,is,ic)  = median(ae(fin));
            end
            if sum(fin) > 1
                sd_all(iL,is,ic) = std(ae(fin));
            end

            fprintf('     %9s  %5s   %11.1f%%   %8.2f  %8.2f     %8.2f   %8.3f\n', ...
                    COLHDR{ic}, SNR_LABEL{is}, res_all(iL,is,ic), ...
                    bias_all(iL,is,ic), sd_all(iL,is,ic), ...
                    med_all(iL,is,ic), spur_all(iL,is,ic));
        end

        if ceil_n(iL,ic) == ntrue
            ctxt = 'resolvable';
        else
            ctxt = sprintf('%d of %d -- NOT resolvable', ceil_n(iL,ic), ntrue);
        end
        % %12s, not %11s: the data rows above carry a literal %% after their
        % %11.1f, so the ceiling row needs one more column to stay aligned.
        fprintf('     %9s  truth   %12s   %8.2f                          %8.3f   <- ceiling, %s\n', ...
                COLHDR{ic}, '--', ceil_err(iL,ic), ...
                max(ceil_n(iL,ic)-ntrue, 0), ctxt);
    end
end
fprintf('\n');
fprintf('   Reading this table: a low "correct count" next to a ceiling that says NOT\n');
fprintf('   resolvable is the angular order failing, not the method. A low count next to\n');
fprintf('   a resolvable ceiling is the method or the noise, and the SNR column says\n');
fprintf('   which -- if it is still low at SNR inf, noise was never the problem.\n\n');

% CHECK. The noise-free arm must recover the true fibre count wherever the
% truth itself resolves. With no noise the only things left between the fit and
% the truth are the band limit and the estimated kernel, so a failure here makes
% nothing at a finite SNR interpretable. Conditions whose ceiling is already
% below the true count are *skipped*, not failed: there the truth does not
% resolve either, and demanding that the fit do better than the truth would be
% the wrong test. (Both manuscript conditions resolve at every Lmax; a 30 degree
% crossing, were one added, would not -- see notebooks/README.md.)
%
% Note also what the noise-free block IS: at sigma = 0 every realisation of a
% condition is the same signal vector, and SMI.fit trains one regression per
% call, so all NREP voxels of that block come back bit identical. That makes its
% std exactly zero -- checked below -- and it means the SNR inf block costs a
% full fit to compute one distinct answer NREP times.
i_inf = find(isinf(SNR_LIST), 1);
if ~isempty(i_inf)
    ok_nf = true; n_tested = 0; ok_det = true;
    for iL = 1:numel(LMAX_LIST)
        for ic = 1:NCOND
            if ceil_n(iL,ic) == numel(axes_gt{ic})
                ok_nf = ok_nf && (res_all(iL,i_inf,ic) > 99);
                n_tested = n_tested + 1;
            end
            ok_det = ok_det && (sd_all(iL,i_inf,ic) == 0);
        end
    end
    fprintf('   CHECK noise-free arm recovers the true count (%d of %d cells resolve)  %s\n', ...
            n_tested, numel(LMAX_LIST)*NCOND, VERDICT{1+ok_nf});
    fprintf('   CHECK noise-free realisations are identical, std exactly 0            %s\n', ...
            VERDICT{1+ok_det});
end

% CHECK. More noise must not help. Compared at the extremes rather than pairwise
% so that Monte Carlo error between adjacent SNRs cannot trip it.
is_lo = SNR_ORD(1); is_hi = SNR_ORD(end);
ok_mono = true;
for iL = 1:numel(LMAX_LIST)
    for ic = 1:NCOND
        ok_mono = ok_mono && (bias_all(iL,is_hi,ic) <= bias_all(iL,is_lo,ic));
    end
end
fprintf('   CHECK mean angular error at SNR %s is <= that at SNR %s   %s\n\n', ...
        SNR_LABEL{is_hi}, SNR_LABEL{is_lo}, VERDICT{1+ok_mono});

%% Step 8 -- export the fODFs so MRtrix can check them
% Everything above is scored by this file. That is fine for reading, but it is
% not independent: the same code wrote the fODF and found its peaks. The point
% of this step is to hand the fODFs to something that has never seen this
% package.
%
% The fit ran at |CS_phase = 0|, where SMI's spherical harmonic basis *is*
% MRtrix's, so the coefficients need no conversion -- they are written straight
% out as MRtrix SH images, one per Lmax.

% The exported image holds the WHOLE sweep -- every SNR block, in the order
% Step 5 laid them out -- so a single sh2peaks run covers all noise levels and
% the key below says which voxel was which.
edir = fullfile(pkgdir, 'export');
if ~exist(edir, 'dir'), mkdir(edir); end
GRID_ALL = C.pick_grid(NVOX);
MR = mrtrix_io();
fprintf('Step 8: writing MRtrix SH images to %s\n', edir);
for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);
    fn = fullfile(edir, sprintf('smifod_lmax%d', Lf));
    MR.write(fn, reshape(fits{iL}.sh, [GRID_ALL size(fits{iL}.sh,2)]));
    fprintf('   smifod_lmax%d.mih  [%s x %d]\n', Lf, mat2str(GRID_ALL), ...
            size(fits{iL}.sh,2));
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
fid = fopen(fullfile(edir,'voxel_key.txt'), 'w');
fprintf(fid, '%% voxel  SNR  crossing_deg  axis1_x axis1_y axis1_z  axis2_x axis2_y axis2_z\n');
fprintf(fid, '%% axis2 is 0 0 0 for the single fibre condition. SNR Inf is the noise-free arm.\n');
fprintf(fid, '%% Voxel order is column-major over the %s grid, matching the .mih images:\n', ...
        mat2str(GRID_ALL));
fprintf(fid, '%% %d contiguous blocks of %d voxels, one per SNR, in the order %s.\n', ...
        NSNR, NVOX_SNR, strjoin(SNR_LABEL, ', '));
fprintf(fid, '%d %g %d %.9f %.9f %.9f %.9f %.9f %.9f\n', key');
fclose(fid);
fprintf('   voxel_key.txt     SNR and true fibre axes per voxel, for matching peaks back\n\n');

fprintf('   To check these against MRtrix3, from %s:\n', edir);
fprintf('     sh2peaks smifod_lmax6.mih peaks_lmax6.mih -num 4\n');
fprintf('     mrinfo   smifod_lmax6.mih\n');
fprintf('     mrconvert smifod_lmax6.mih smifod_lmax6.mif   # if you prefer a single file\n');
fprintf('     mrview   smifod_lmax6.mih -odf.load_sh smifod_lmax6.mih\n');
fprintf('   sh2peaks writes 3 components per peak; compare against voxel_key.txt.\n');
fprintf('   Nothing in this walkthrough runs those commands -- that is the point.\n\n');

%% Figures
% The six manuscript figures. Set |MAKE_FIGURES = false| to skip them.
%
% All glyphs use |SMI_response_helpers|, the renderer MRtrix's |shview| logic
% implies: *radius is |amplitude|, colour is the signed amplitude*, so negative
% lobes show as a colour change instead of being folded silently into the
% surface. That matters here because the band-limited truth really is negative
% over much of the sphere (Step 2).
%
% Every 3D panel opens in an *isometric* view, so all subfigures are readable
% without touching the camera:
%
%  Fig 1  the two ground truth fibre configurations
%  Fig 2  a montage of zonal harmonic response glyphs, b = 0 to 3
%  Fig 3  the spherical signal: *b down the rows*, SNR across the columns,
%         increasing left to right and ending at the ground truth
%  Fig 4  reconstructed fODFs, *single fibre*: Lmax down the rows, the ground
%         truth in the first column and then one column per SNR
%  Fig 5  reconstructed fODFs, *60 degree crossing*: the same layout
%  Fig 6  bias, spread and spurious peaks against SNR -- the summary plots
%
% Figures 3 to 5 all run SNR left to right in *increasing* order and put the
% noiseless reference in the *last* column of Figure 3 and the *first* column of
% Figures 4 and 5, next to the axis it belongs against in each case: the signal
% degrades away from the truth, and the reconstructions are read against it.

MAKE_FIGURES = true;
ISO_VIEW     = [1 1 1];      % isometric camera: equal foreshortening on x, y, z
GLYPH_N      = 121;          % mesh resolution for every glyph

if MAKE_FIGURES
    fprintf('Figures: rendering ...\n');
    [THg, PHg, dirs_g] = RH.grid(GLYPH_N, GLYPH_N);

    % ================================================================
    % FIGURE 1 -- the two ground truth fibre configurations, isometric
    % ================================================================
    % Both panels open already rotated to an isometric view, so the 60 degree
    % crossing reads as a crossing without anyone having to drag the camera.
    figure('Name', 'Fig 1  ground truth fibre geometry');
    for ic = 1:NCOND
        subplot(1, NCOND, ic);
        [X, Y, Z, Cc] = RH.sh_glyph(plm_gt(ic,:), LMAX_GT, C.CS_PHASE, ...
                                    GLYPH_N, GLYPH_N, 1);
        surf(X, Y, Z, Cc); shading interp;
        axis equal off vis3d; view(ISO_VIEW);
        camlight headlight; lighting gouraud;
        title(sprintf('%s   (kappa = %g, Lmax %d)', COLHDR{ic}, C.KAPPA, LMAX_GT));
    end

    % ================================================================
    % FIGURE 2 -- montage of zonal harmonic response glyphs, b = 0 to 3
    % ================================================================
    % One glyph per shell, drawn by the same renderer as Figure 1 so a response
    % and an fODF can be compared without a convention change in between. This
    % is the panel the CSD and MSMT-CSD responses will be added to: their
    % dwi2response output is zonal coefficients in exactly this form, so each
    % becomes another row here.
    %
    % Each glyph is normalised to its OWN maximum, so the montage compares
    % SHAPE across shells. The amplitude falls steeply with b -- that is the
    % r_0 column of the table in Step 3, not something to read off the glyphs.
    figure('Name', 'Fig 2  response glyph per shell');
    nsh = numel(b_shell);
    for i = 1:nsh
        subplot(1, nsh, i);
        r  = RH.zh(K, b_shell(i), LMAX_GT, C.D_FW);
        pk = max(abs(RH.profile(r, linspace(0, pi, 361))));
        if pk <= 0, pk = 1; end
        [X, Y, Z, Cc] = RH.zh_glyph(r, GLYPH_N, GLYPH_N, 1/pk);
        surf(X, Y, Z, Cc); shading interp;
        axis equal off vis3d; view(ISO_VIEW);
        camlight headlight; lighting gouraud;
        title(sprintf('b = %.0f', b_shell(i)));
    end

    % ================================================================
    % FIGURE 3 -- the spherical signal: b down the rows, SNR across
    % ================================================================
    % What the forward convolution actually produces, as a surface on the
    % sphere rather than as numbers: radius is S(u)/S0.
    %
    % *b is the vertical axis*, increasing downwards, so a column is one noise
    % level read down the shells and a row is one shell read across the noise
    % levels. SNR increases left to right and the last column is the ground
    % truth: the model evaluated directly on a dense grid, with no noise and no
    % projection.
    %
    % The SNR columns cannot be drawn that way -- only 90 directions were
    % measured -- so each is the spherical harmonic fit of ONE representative
    % realisation at that SNR, at the same Lmax the ground truth uses. That is
    % the same projection SMI does internally, so those surfaces are what the
    % fit actually sees. The |inf| column and the ground truth column are
    % therefore *different computations of the same object*, and they should be
    % indistinguishable; that they are is a check you can make by eye.
    %
    % Every panel shares one radial scale, so the fall in signal with b is
    % visible rather than normalised away.
    figure('Name', 'Fig 3  spherical signal: b down the rows, SNR across the columns');
    ic_show = NCOND;                                  % the 60 degree crossing
    dwsh    = 2:nsh;                                  % skip the b~0 shell
    nrow3   = numel(dwsh);
    ncol3   = NSNR + 1;                               % every SNR, then the truth
    Yg      = SMI.get_even_SH(dirs_g, LMAX_GT, C.CS_PHASE);

    amps = cell(nrow3, ncol3);
    for j = 1:nrow3
        i = dwsh(j);
        m = (shell_id(:)' == i);
        Ym = SMI.get_even_SH(bvecs(m,:), LMAX_GT, C.CS_PHASE);
        for k = 1:NSNR
            is = SNR_ORD(k);
            v  = find(cond_id == ic_show & snr_id == is, 1);
            cm = Ym \ S_noisy(v, m)';                 % least squares SH fit
            amps{j,k} = reshape(Yg*cm, size(THg));
        end
        bq = b_shell(i)*ones(1, size(dirs_g,1));
        s_clean_g = H.signal(plm_gt(ic_show,:), [K 1 1], bq, ones(size(bq)), ...
                             zeros(size(bq)), dirs_g, LMAX_GT, C.CS_PHASE, C.D_FW);
        amps{j,ncol3} = reshape(s_clean_g, size(THg));
    end
    smax = 0;
    for k = 1:numel(amps), smax = max(smax, max(abs(amps{k}(:)))); end

    for j = 1:nrow3
        for k = 1:ncol3
            subplot(nrow3, ncol3, (j-1)*ncol3 + k);
            [X, Y, Z, Cc] = RH.glyph(amps{j,k}, THg, PHg, 1/smax);
            surf(X, Y, Z, Cc); shading interp;
            axis equal off vis3d; view(ISO_VIEW);
            camlight headlight; lighting gouraud;
            if k <= NSNR
                title(sprintf('b = %.0f, SNR %s', b_shell(dwsh(j)), ...
                              SNR_LABEL{SNR_ORD(k)}));
            else
                title(sprintf('b = %.0f, ground truth', b_shell(dwsh(j))));
            end
        end
    end

    % ================================================================
    % FIGURES 4 and 5 -- reconstructed fODFs, one figure per condition
    % ================================================================
    % Same renderer and same camera as Figure 1, so a reconstruction can be put
    % next to the truth and compared directly -- which is now literal: the first
    % column of each figure IS the Figure 1 glyph for that condition, truncated
    % to the row's Lmax.
    %
    % Truncating it per row rather than always drawing it at |LMAX_GT| is the
    % honest comparison, and it is the same object Step 7 scores the ceiling on:
    % at Lmax 4 and 6 the truth carries detail the fit cannot represent, and the
    % first column shows exactly how much. At Lmax 8, where |LMAX_GT| is also 8,
    % it is Figure 1's panel unchanged.
    %
    % Columns after the first are one realisation per SNR, increasing left to
    % right, so a row reads as "the truth, then what came back out of the noise".
    for ic = 1:NCOND
        figure('Name', sprintf('Fig %d  reconstructed fODFs, %s', 3+ic, COLHDR{ic}));
        for iL = 1:numel(LMAX_LIST)
            Lf = LMAX_LIST(iL);
            nc = (Lf/2+1)*(Lf+1);
            Lv = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
            scl = sqrt((2*Lv(2:end)'+1)/(4*pi));

            % column 1: the ground truth, truncated to this Lmax
            subplot(numel(LMAX_LIST), NSNR+1, (iL-1)*(NSNR+1) + 1);
            [X, Y, Z, Cc] = RH.sh_glyph(plm_gt(ic, 1:nc-1), Lf, C.CS_PHASE, ...
                                        GLYPH_N, GLYPH_N, 1);
            surf(X, Y, Z, Cc); shading interp;
            axis equal off vis3d; view(ISO_VIEW);
            camlight headlight; lighting gouraud;
            title(sprintf('truth, Lmax %d', Lf));

            % then one column per SNR, increasing left to right
            for k = 1:NSNR
                is = SNR_ORD(k);
                subplot(numel(LMAX_LIST), NSNR+1, (iL-1)*(NSNR+1) + 1 + k);
                v = find(cond_id == ic & snr_id == is, 1);
                % back to the p_00 = 1 convention the glyph renderer expects
                plm_v = fits{iL}.sh(v, 2:nc) ./ scl;
                [X, Y, Z, Cc] = RH.sh_glyph(plm_v, Lf, C.CS_PHASE, ...
                                            GLYPH_N, GLYPH_N, 1);
                surf(X, Y, Z, Cc); shading interp;
                axis equal off vis3d; view(ISO_VIEW);
                camlight headlight; lighting gouraud;
                title(sprintf('Lmax %d, SNR %s', Lf, SNR_LABEL{is}));
            end
        end
    end

    % ================================================================
    % FIGURE 6 -- bias, spread and spurious peaks against SNR
    % ================================================================
    % *Preliminary.* These are the summary plots the sweep exists to produce,
    % and they are drawn from the same arrays Step 7 printed, so the figure and
    % the table cannot disagree.
    %
    % Three things worth knowing before reading them:
    %
    % * *"Bias" here is the mean angular error*, which is the mean of a
    %   non-negative quantity and so does not vanish for a perfect estimator.
    %   The floor is set by the direction grid (about 1.5 degrees at 1500
    %   directions) and by the band limit. The right way to read it is against
    %   the noise-free point of the same curve, not against zero.
    % * *The x axis is the SNR list in increasing order, evenly spaced*, not to
    %   scale. |inf| cannot go on a numeric axis, and equal spacing keeps the
    %   low-SNR end -- where everything happens -- from being squeezed.
    % * *Spurious peak counts and angular error fail differently.* A method can
    %   hold its angular error and start inventing peaks, or lose the peak
    %   entirely and report a confident wrong direction. Both panels are needed.
    figure('Name', 'Fig 6  bias, spread and spurious peaks vs SNR');
    mks    = {'-o', '-s', '-^', '-d'};
    mets   = {bias_all, sd_all, spur_all};
    mnames = {'bias: mean angular error (deg)', ...
              'std of angular error (deg)', ...
              'spurious peaks per voxel'};
    lgd = cell(1, numel(LMAX_LIST));
    for iL = 1:numel(LMAX_LIST)
        lgd{iL} = sprintf('Lmax %d', LMAX_LIST(iL));
    end
    for ic = 1:NCOND
        for im = 1:numel(mets)
            subplot(NCOND, numel(mets), (ic-1)*numel(mets) + im);
            hold on;
            M = mets{im};
            for iL = 1:numel(LMAX_LIST)
                y = M(iL, SNR_ORD, ic);
                plot(1:NSNR, y(:)', mks{mod(iL-1, numel(mks))+1});
            end
            set(gca, 'XTick', 1:NSNR, 'XTickLabel', SNR_LABEL(SNR_ORD));
            xlim([0.5 NSNR+0.5]);
            xlabel('SNR'); ylabel(mnames{im});
            title(sprintf('%s -- %s', COLHDR{ic}, mnames{im}));
            grid on;
            if ic == 1 && im == 1, legend(lgd); end
            hold off;
        end
    end

    fprintf('   6 figures drawn; every 3D panel opens in an isometric view.\n');
    fprintf('   Radius is |amplitude| and colour is the SIGNED amplitude, so the\n');
    fprintf('   negative lobes of a band-limited fODF show as a colour change\n');
    fprintf('   rather than being folded into the surface.\n\n');
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
% *What to change first if you want to probe it.* |NREP|, |SNR_LIST| and
% |LMAX_LIST| at the top of this file -- |NREP| is the one that sets the
% runtime. Then |mc_config.m|: |PROTOCOL| names the acquisition, |KAPPA| sets
% fibre dispersion, |K_WM| the tissue, |ANGLES| the crossings. Changing them
% there changes them for the pipeline too, which is the point of that file
% existing.

fprintf('=== walkthrough complete ===\n');
fprintf('%d realisations per condition per SNR, %d SNR levels, Lmax %s\n\n', ...
        NREP, NSNR, mat2str(LMAX_LIST));

sumnames = {'correct fibre count (%)', ...
            'bias: mean angular error (deg)', ...
            'std of angular error (deg)', ...
            'spurious peaks per voxel'};
sums = {res_all, bias_all, sd_all, spur_all};
for im = 1:numel(sums)
    fprintf('%s\n', sumnames{im});
    for ic = 1:NCOND
        fprintf('  %s\n', COLHDR{ic});
        fprintf('        SNR ');
        for k = 1:NSNR, fprintf('%10s', SNR_LABEL{SNR_ORD(k)}); end
        fprintf('\n');
        for iL = 1:numel(LMAX_LIST)
            fprintf('     Lmax %d ', LMAX_LIST(iL));
            for k = 1:NSNR
                fprintf('%10.3f', sums{im}(iL, SNR_ORD(k), ic));
            end
            fprintf('\n');
        end
    end
    fprintf('\n');
end
