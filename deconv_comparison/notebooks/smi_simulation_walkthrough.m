%% Simulating SMI: a step-by-step walkthrough
% This is the SMI arm of the Monte Carlo behind
% |Reports/REPORT_SMI_deconvolution_MonteCarlo.md|, taken apart so that every
% step can be inspected and checked rather than taken on trust.
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
%  Step 2  the ground truth fibre geometry  -> Watson fODFs at 0/15/45/60 deg
%  Step 3  the kernel, and it as a response -> K_l(b), zonal harmonics
%  Step 4  forward convolution              -> noise-free signal
%  Step 5  Rician noise                     -> the measured data
%  Step 6  SMI.fit at Lmax 4, 6 and 8       -> kernel and fODF recovery
%  Step 7  peaks and angular error          -> did it work, and at which Lmax?
%  Step 8  export the fODFs for MRtrix      -> check them outside this file
%  Step 9  how to scale this up
%
% *Runtime is about six minutes*, almost all of it Step 6. |SMI.fit| is
% dominated by training its regression rather than by the voxel count, so each
% Lmax costs roughly the same two minutes whether there are 100 voxels or
% 10,000 -- which is also why the full campaign is affordable.

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

NREP      = 25;                 % realisations per condition. The report used 10000.
SNR       = 50;                 % 1/sigma. The report sweeps 5, 10, 20, 30, 50, inf.
LMAX_LIST = [4 6 8];            % angular orders to fit at -- see the next section
PROTOCOL  = 'hcp_real_3shell.txt';

C = mc_config();                % conditions, kernel, dispersion, seeds
H = fODF_modulation_helpers();
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
fprintf('\n=== SMI simulation walkthrough ===\n');
fprintf('protocol %s, NREP %d per condition, SNR %g, CS_phase %d\n', ...
        PROTOCOL, NREP, SNR, C.CS_PHASE);
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

proto = fileread(fullfile(pkgdir, 'protocol', PROTOCOL));
cols  = textscan(proto, '%f %f %f %f', 'CommentStyle', '%');
bvals = cols{1}(:)';
bvecs = [cols{2} cols{3} cols{4}];
Ndwi  = numel(bvals);

fprintf('Step 1: %d volumes, %d distinct b values\n', Ndwi, numel(unique(bvals)));

% The directions are only unit to the precision the .bvec file was written at.
% That sounds like a nuisance and is not: anything that treats g(3) as cos(theta)
% without normalising inherits the error, and at Lmax 8 that is enough to break
% the Step 3 identity by 5e-7 -- six orders of magnitude worse than it should be.
% So normalise here, once, and report what was corrected.
e_raw = max(abs(sqrt(sum(bvecs.^2, 2)) - 1));
bvecs = bvecs ./ repmat(sqrt(sum(bvecs.^2, 2)), 1, 3);
fprintf('   .bvec unit-norm error as supplied: %.2e, normalised to %.2e\n', ...
        e_raw, max(abs(sqrt(sum(bvecs.^2, 2)) - 1)));

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

% CHECK 1. Every direction must be a unit vector. Now trivially true, but the
% number printed above is the one that matters.
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
% 15, 45 and 60 degrees. Each population is a *Watson* distribution with
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
fprintf('        b   ');  fprintf('%10s', 'single', '15 deg', '45 deg', '60 deg'); fprintf('\n');
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

%% Step 5 -- Rician noise
% Complex Gaussian noise is added to a real signal and the magnitude taken,
% which is exactly Rician:
%
%  S_noisy = sqrt( (S + sigma*n1)^2 + (sigma*n2)^2 ),   n1, n2 ~ N(0,1)
%
% |sigma = 1/SNR|, with SNR defined against the b = 0 signal. The seed comes
% from |mc_config.m|, so this walkthrough and the pipeline draw the same noise.

sigma   = 1/SNR;
NVOX    = NCOND*NREP;
GRID    = C.pick_grid(NVOX);
cond_id = repelem((1:NCOND)', NREP, 1);
S_rep   = S_clean(cond_id, :);

rand('seed', C.SEED_MC); randn('seed', C.SEED_MC);
S_noisy = sqrt((S_rep + sigma*randn(size(S_rep))).^2 + ...
               (        sigma*randn(size(S_rep))).^2);
dwi = reshape(S_noisy, [GRID Ndwi]);

fprintf('Step 5: %d conditions x %d reps = %d voxels, grid %s\n', ...
        NCOND, NREP, NVOX, mat2str(GRID));
fprintf('   sigma = 1/SNR = %.4f\n', sigma);
fprintf('   the truth''s band-limiting error is %.0fx smaller than sigma, so the\n', ...
        sigma/e_sh);
fprintf('   band limit is not what limits any result below\n');

% CHECK 1. Recover sigma from the simulated data rather than trusting the value
% typed in. This is an independent read of the noise level.
resid = S_noisy - S_rep;
s_hat = std(resid(:));
rel   = abs(s_hat - sigma)/sigma;
fprintf('   CHECK recovered sigma       %.5f vs %.5f, %.1f%% off   %s\n', ...
        s_hat, sigma, 100*rel, VERDICT{1+(rel < 0.10)});

% CHECK 2. A magnitude is strictly positive, which is why SMI is told
% NoiseBias = 'Rician' in Step 6.
mn = min(S_noisy(:));
fprintf('   CHECK signal strictly positive   min = %.4f            %s\n', ...
        mn, VERDICT{1+(mn > 0)});
hib = (shell_id(:)' == numel(b_shell));
fprintf('   Rician floor at b = %.0f: mean %.4f against a noise-free %.4f (upward bias)\n\n', ...
        b_shell(end), mean(mean(S_noisy(:,hib))), mean(mean(S_rep(:,hib))));

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
% The whole fit is repeated for each Lmax, so everything downstream can be read
% as a function of angular order.

fits = cell(1, numel(LMAX_LIST));
fprintf('Step 6: SMI.fit on %d voxels, once per Lmax\n', NVOX);
for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);

    options = struct();
    options.b     = bvals;
    options.dirs  = bvecs;
    options.sigma = sigma*ones(GRID);
    options.mask  = true(GRID);
    options.compartments  = {'IAS','EAS','FW'};
    options.NoiseBias     = 'Rician';
    options.Lmax          = [0 Lf Lf Lf];
    options.CS_phase      = C.CS_PHASE;
    options.D_FW          = C.D_FW;
    options.flag_fit_fODF = 1;
    options.fODF_regularization = struct('flag_nonneg', 1, 'lambda_tikhonov', 0.3);

    t0  = tic;
    out = SMI.fit(dwi, options);
    el  = toc(t0);

    Lv  = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
    plm = reshape(out.plm, [NVOX numel(Lv)-1]);
    sh  = [ones(NVOX,1) plm] .* repmat(sqrt((2*Lv+1)/(4*pi))', NVOX, 1);
    nonfinite = sum(~isfinite(sh(:)));
    sh(~isfinite(sh)) = 0;

    fits{iL} = struct('Lmax', Lf, 'sh', sh, 'nonfinite', nonfinite, ...
                      'kernel', reshape(out.kernel, [NVOX size(out.kernel,4)]), ...
                      'converged', out.fODF_regularization.flag_converged(:), ...
                      'seconds', el);
    fprintf('   Lmax %d: %2d coefficients, %5.1f s  (lambda_nonneg = %g, lambda_tikhonov = %g)\n', ...
            Lf, numel(Lv), el, out.fODF_regularization.lambda_nonneg, ...
            out.fODF_regularization.lambda_tikhonov);
end

% The kernel is estimated from rotational invariants, which do not depend on
% the fODF's Lmax in the same way -- so this table should be nearly flat across
% Lmax. If it is not, the kernel fit is being disturbed by the fODF fit.
knm = {'f','Da','Depar','Deperp','fw'};
fprintf('   kernel recovery, median over %d voxels (true value in brackets):\n', NVOX);
fprintf('        Lmax  ');
for j = 1:5, fprintf('%16s', sprintf('%s [%.2f]', knm{j}, K(j))); end
fprintf('\n');
for iL = 1:numel(LMAX_LIST)
    fprintf('        %4d  ', LMAX_LIST(iL));
    for j = 1:5
        col = fits{iL}.kernel(:,j); col = col(isfinite(col));
        fprintf('%16.3f', median(col));
    end
    fprintf('\n');
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
fprintf('   CHECK fODF all finite at every Lmax                    %s\n', VERDICT{1+all_fin});

% CHECK. The p_00 = 1 convention must survive the fit at every Lmax.
e_l0 = 0;
for iL = 1:numel(LMAX_LIST)
    e_l0 = max(e_l0, max(abs(fits{iL}.sh(:,1) - 1/sqrt(4*pi))));
end
fprintf('   CHECK p_00 == 1 convention held   max|err| = %.2e   %s\n\n', ...
        e_l0, VERDICT{1+(e_l0 < 1e-12)});

%% Step 7 -- peaks and angular error, per Lmax
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
% method failed" from "this angular order cannot represent the answer".

PEAK_NBR = 12;      % degrees, neighbourhood for the local-maximum test
PEAK_REL = 0.30;    % keep peaks at >= 30% of the largest anisotropic amplitude
de  = H.dirs(1500);
nbr = (de*de') > cosd(PEAK_NBR);
nbr(1:size(de,1)+1:end) = true;               % include self, so >= means "is the max"

fprintf('Step 7: peaks (%d directions, %d deg neighbourhood, %.0f%% threshold)\n', ...
        size(de,1), PEAK_NBR, 100*PEAK_REL);
fprintf('   The grid spacing sets a floor on angular error of roughly %.1f deg;\n', ...
        sqrt(4*pi/size(de,1))*180/pi/2);
fprintf('   the ceiling column shows that floor directly.\n\n');
fprintf('   Lmax  condition   correct count   median error   spurious   ceiling\n');

res_all = zeros(numel(LMAX_LIST), NCOND);
for iL = 1:numel(LMAX_LIST)
    Lf   = LMAX_LIST(iL);
    Ye   = SMI.get_even_SH(de, Lf, C.CS_PHASE);
    nc   = size(Ye,2);
    Afit = fits{iL}.sh * Ye' - 1/(4*pi);
    Agt  = sh_gt(:,1:nc) * Ye' - 1/(4*pi);    % truth truncated to THIS Lmax

    for ic = 1:NCOND
        rows   = find(cond_id == ic);
        ntrue  = numel(axes_gt{ic});
        nfound = zeros(numel(rows),1);
        aerr   = nan(numel(rows),1);
        for r = 1:numel(rows)
            a  = Afit(rows(r),:)';
            lm = find(a > 0 & a >= max(a .* nbr, [], 1)');
            if isempty(lm), continue; end
            lm = lm(a(lm) >= PEAK_REL*max(a(lm)));
            [~, o] = sort(a(lm), 'descend'); lm = lm(o);
            P = de(lm,:);
            sel = true(size(P,1),1);
            for i = 1:size(P,1)
                if ~sel(i), continue; end
                dup = abs(P*P(i,:)') > cosd(PEAK_NBR); dup(i) = false; sel(dup) = false;
            end
            P = P(sel,:);
            nfound(r) = size(P,1);
            d = zeros(1,ntrue);
            for k = 1:ntrue, d(k) = acosd(min(abs(P(1,:)*axes_gt{ic}{k}'), 1)); end
            aerr(r) = min(d);
        end

        % the ceiling for this Lmax and condition
        a  = Agt(ic,:)';
        lm = find(a > 0 & a >= max(a .* nbr, [], 1)');
        lm = lm(a(lm) >= PEAK_REL*max(a(lm)));
        [~, o] = sort(a(lm), 'descend'); lm = lm(o);
        P = de(lm,:);
        sel = true(size(P,1),1);
        for i = 1:size(P,1)
            if ~sel(i), continue; end
            dup = abs(P*P(i,:)') > cosd(PEAK_NBR); dup(i) = false; sel(dup) = false;
        end
        nceil = size(P(sel,:),1);

        res_all(iL,ic) = 100*mean(nfound == ntrue);
        if nceil == ntrue, ctxt = 'resolvable';
        else,              ctxt = sprintf('%d of %d -- NOT resolvable', nceil, ntrue);
        end
        fprintf('   %4d  %6d deg   %11.1f%%   %9.2f deg   %8.3f   %s\n', ...
                Lf, C.ANGLES(ic), res_all(iL,ic), median(aerr(isfinite(aerr))), ...
                mean(max(nfound - ntrue, 0)), ctxt);
    end
end
fprintf('\n');
fprintf('   Reading this table: a low "correct count" next to a ceiling that says NOT\n');
fprintf('   resolvable is the angular order failing, not the method. A low count next to\n');
fprintf('   a resolvable ceiling is the method or the noise.\n\n');

%% Step 8 -- export the fODFs so MRtrix can check them
% Everything above is scored by this file. That is fine for reading, but it is
% not independent: the same code wrote the fODF and found its peaks. The point
% of this step is to hand the fODFs to something that has never seen this
% package.
%
% The fit ran at |CS_phase = 0|, where SMI's spherical harmonic basis *is*
% MRtrix's, so the coefficients need no conversion -- they are written straight
% out as MRtrix SH images, one per Lmax.

edir = fullfile(pkgdir, 'export');
if ~exist(edir, 'dir'), mkdir(edir); end
MR = mrtrix_io();
fprintf('Step 8: writing MRtrix SH images to %s\n', edir);
for iL = 1:numel(LMAX_LIST)
    Lf = LMAX_LIST(iL);
    fn = fullfile(edir, sprintf('smifod_lmax%d', Lf));
    MR.write(fn, reshape(fits{iL}.sh, [GRID size(fits{iL}.sh,2)]));
    fprintf('   smifod_lmax%d.mih  [%s x %d]\n', Lf, mat2str(GRID), size(fits{iL}.sh,2));
end

% A voxel-order key, so a peak found in MRtrix can be matched to the condition
% that produced it. Column 1 is the linear voxel index in MRtrix's order,
% column 2 the condition, columns 3-5 and 6-8 the true fibre axes.
key = zeros(NVOX, 8);
for v = 1:NVOX
    ic = cond_id(v);
    key(v,1) = v; key(v,2) = C.ANGLES(ic);
    key(v,3:5) = axes_gt{ic}{1};
    if numel(axes_gt{ic}) == 2, key(v,6:8) = axes_gt{ic}{2}; end
end
fid = fopen(fullfile(edir,'voxel_key.txt'), 'w');
fprintf(fid, '%% voxel  crossing_deg  axis1_x axis1_y axis1_z  axis2_x axis2_y axis2_z\n');
fprintf(fid, '%% axis2 is 0 0 0 for the single fibre condition.\n');
fprintf(fid, '%% Voxel order is column-major over the %s grid, matching the .mih images.\n', ...
        mat2str(GRID));
fprintf(fid, '%d %d %.9f %.9f %.9f %.9f %.9f %.9f\n', key');
fclose(fid);
fprintf('   voxel_key.txt     true fibre axes per voxel, for matching peaks back\n\n');

fprintf('   To check these against MRtrix3, from %s:\n', edir);
fprintf('     sh2peaks smifod_lmax6.mih peaks_lmax6.mih -num 4\n');
fprintf('     mrinfo   smifod_lmax6.mih\n');
fprintf('     mrconvert smifod_lmax6.mih smifod_lmax6.mif   # if you prefer a single file\n');
fprintf('     mrview   smifod_lmax6.mih -odf.load_sh smifod_lmax6.mih\n');
fprintf('   sh2peaks writes 3 components per peak; compare against voxel_key.txt.\n');
fprintf('   Nothing in this walkthrough runs those commands -- that is the point.\n\n');

%% Step 9 -- from here to a full campaign
% Everything above is one SNR at |NREP| realisations. The published campaign is
% the same forward model at |NREP = 10000| across six noise levels, driven by
% |gen_montecarlo.m| rather than by hand:
%
%  cd deconv_comparison
%  octave --eval "run('oct_path.m'); gen_montecarlo(50, 10000, 'snr50')"
%  ./run_mrtrix.sh fit snr50            % CSD, MSMT-CSD and every peak, in MRtrix
%  python3 tables.py nf snr50:50 ...    % the tables in Reports/deconv_tables.md
%
% *Note that the pipeline still reads the older synthetic protocol.* This
% walkthrough runs on the real HCP scheme; |gen_montecarlo.m| has not been
% switched over, so the published tables and this file are not yet measuring
% the same acquisition. Switching it is a one-line change in |gen_montecarlo.m|
% plus a full re-run.
%
% *What this walkthrough leaves out on purpose.* The comparison arms. CSD and
% MSMT-CSD are MRtrix3 binaries driven by |run_mrtrix.sh|. This file is the SMI
% arm alone, taken apart.
%
% *What to change first if you want to probe it.* |NREP|, |SNR|, |LMAX_LIST|
% and |PROTOCOL| at the top of this file. Then |mc_config.m|: |KAPPA| sets
% fibre dispersion, |K_WM| the tissue, |ANGLES| the crossings.

fprintf('=== walkthrough complete ===\n');
fprintf('correct fibre count, %% of %d realisations at SNR %g:\n', NREP, SNR);
fprintf('        ');  fprintf('%10s', 'single', '15 deg', '45 deg', '60 deg'); fprintf('\n');
for iL = 1:numel(LMAX_LIST)
    fprintf('  Lmax %d', LMAX_LIST(iL)); fprintf('%10.1f', res_all(iL,:)); fprintf('\n');
end
fprintf('\n');
