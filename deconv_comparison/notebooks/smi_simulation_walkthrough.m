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
% |helpers/fODF_modulation_helpers.m|, the fit is the real |SMI.fit|, and every
% constant comes from |mc_config.m| -- the same file |gen_montecarlo.m| reads.
% If a number here disagrees with the pipeline, the pipeline is what is wrong.
%
% *What you should end up believing.* Not "the code ran". Each step ends with
% one or more *CHECK* lines that compare its output against something computed
% a different way, and print |ok| or |** FAILED **| next to the number. If you
% read only the CHECK lines, you have audited the simulation.
%
%  Step 1  the acquisition protocol         -> shells and directions
%  Step 2  the ground truth fibre geometry  -> Watson fODFs at 0/15/45/60 deg
%  Step 3  the kernel, and it as a response -> K_l(b), zonal harmonics
%  Step 4  forward convolution              -> noise-free signal
%  Step 5  Rician noise                     -> the measured data
%  Step 6  SMI.fit                          -> kernel and fODF recovery
%  Step 7  peaks and angular error          -> did it work?
%  Step 8  cross-check the peaks against MRtrix (skipped if not installed)
%  Step 9  how to scale this up to the published run
%
% Runtime is a couple of minutes at the default size.

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

NREP = 50;      % realisations per condition. 50 is quick; the report used 10000.
SNR  = 50;      % 1/sigma. The report sweeps 5, 10, 20, 30, 50 and noise free.

C = mc_config();              % the single definition of the experiment
H = fODF_modulation_helpers();
VERDICT = {'** FAILED **', 'ok'};             % VERDICT{1+condition}

fprintf('\n=== SMI simulation walkthrough ===\n');
fprintf('NREP %d per condition, SNR %g, Lmax fit %d, Lmax truth %d, CS_phase %d\n\n', ...
        NREP, SNR, C.LMAX_FIT, C.LMAX_GT, C.CS_PHASE);

%% Step 1 -- the acquisition protocol
% Three shells at b = 1, 2, 3 ms/um^2 with 90 directions each, plus 18 b = 0:
% the HCP-like protocol Jeurissen et al. (2014) used to introduce MSMT-CSD.
%
% The directions were generated once by electrostatic repulsion
% (|setup_protocol.py|, the only file in the package that needs dipy) and are
% *tracked as text* at |protocol/hcp_like_3shell.txt|, so this walkthrough
% needs no Python and reads the identical table every other arm reads. b is
% carried in ms/um^2 throughout; b = 1 here is 1000 s/mm^2.

proto = fileread(fullfile(pkgdir, 'protocol', 'hcp_like_3shell.txt'));
cols  = textscan(proto, '%f %f %f %f', 'CommentStyle', '%');
bvals = cols{1}(:)';
bvecs = [cols{2} cols{3} cols{4}];
Ndwi  = numel(bvals);
shells = unique(bvals);          % already a row: bvals is a row

fprintf('Step 1: %d volumes\n', Ndwi);
for s = shells
    fprintf('   b = %g : %3d directions\n', s, sum(bvals == s));
end

% CHECK. The diffusion-weighted directions must be unit vectors, and the b = 0
% rows must carry no direction at all. The forward model assumes both silently,
% so they are worth making loud.
dw    = bvals > 0;
gdw   = bvecs(dw,:);
gb0   = bvecs(~dw,:);
e_nrm = max(abs(sqrt(sum(gdw.^2, 2)) - 1));
e_b0  = max(abs(gb0(:)));
fprintf('   CHECK directions are unit   max| |g|-1 | = %.2e   %s\n', ...
        e_nrm, VERDICT{1+(e_nrm < 1e-12)});
fprintf('   CHECK b=0 rows carry no dir max|g|       = %.2e   %s\n\n', ...
        e_b0, VERDICT{1+(e_b0 == 0)});

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
L_gt  = repelem(0:2:C.LMAX_GT, 2*(0:2:C.LMAX_GT)+1)';
keepF = L_gt <= C.LMAX_FIT;

fodf_gt = zeros(size(dq,1), NCOND);
plm_gt  = zeros(NCOND, numel(L_gt)-1);
sh_gt   = zeros(NCOND, sum(keepF));
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
    p             = H.mixture_plm(f, dq, C.LMAX_GT, C.CS_PHASE);
    plm_gt(ic,:)  = p(:)';
    coef          = [1; p(:)] .* sqrt((2*L_gt+1)/(4*pi));
    sh_gt(ic,:)   = coef(keepF)';

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
Yq    = SMI.get_even_SH(dq, C.LMAX_GT, C.CS_PHASE);
sc_gt = sqrt((2*L_gt' + 1)/(4*pi));
amp   = ([ones(NCOND,1) plm_gt] .* repmat(sc_gt, NCOND, 1)) * Yq';
mass  = mean(amp, 2) * 4*pi;                  % equal-area quadrature
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
% Watson mixture at Lmax rings, exactly as truncating any Fourier series does,
% and the rings go below zero. So the "truth" the non-negativity constraint in
% Step 6 is asked to approach does not satisfy non-negativity either. This is
% why the constraint is a regularizer and not a statement of fact, and why the
% regularization report quotes the ground truth's own negative mass whenever it
% quotes an estimate's.
fprintf('   band-limited truth at Lmax %d -- peak, minimum, and %% of the sphere below 0:\n', ...
        C.LMAX_GT);
for ic = 1:NCOND
    fprintf('     %2d deg : peak %+.4f, min %+.4f, negative over %4.1f%% of directions\n', ...
            C.ANGLES(ic), max(amp(ic,:)), min(amp(ic,:)), 100*mean(amp(ic,:) < 0));
end
fprintf('   The ringing is shallow next to the peak, but it covers a lot of the sphere,\n');
fprintf('   because the fODF is near zero almost everywhere and the rings straddle zero.\n');
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
fprintf('   zonal response r_l, one row per shell -- this IS an MRtrix response file:\n');
fprintf('        b        r_0       r_2       r_4       r_6\n');
r_all = RH.zh(K, shells, C.LMAX_FIT, C.D_FW);       % [Nshell x (Lmax/2+1)]
for i = 1:numel(shells)
    fprintf('     %4g  %9.4f %9.4f %9.4f %9.4f\n', shells(i), r_all(i,:));
end

% CHECK. Build the signal of a delta fODF along z two ways and compare.
%
% Route A, the zonal profile above. Route B, SMI's own forward model fed the
% plm of a delta. In the p_00 = 1 convention a delta along z is exact and needs
% no approximation: f_lm = Y_lm(z), which is zero for every m /= 0 and
% sqrt((2l+1)/4pi) at m = 0, so every p_l0 = 1 and every other p_lm = 0.
Mf = [];
for l = 2:2:C.LMAX_FIT, Mf = [Mf, -l:l]; end
plm_delta = double(Mf(:)' == 0);               % 1 at m = 0, 0 elsewhere

zax   = [0 0 1];
bdw   = bvals(dw);
S_fwd = H.signal(plm_delta, [K 1 1], bdw, ones(1,sum(dw)), zeros(1,sum(dw)), ...
                 gdw, C.LMAX_FIT, C.CS_PHASE, C.D_FW);
theta = acos(max(min(gdw * zax', 1), -1));
S_zon = zeros(sum(dw), 1);
for i = 1:numel(shells)
    m = (bdw == shells(i));
    if any(m), S_zon(m) = RH.profile(r_all(i,:), theta(m)); end
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
                 bvecs, C.LMAX_GT, C.CS_PHASE, C.D_FW);
    S_clean(ic,:) = s(:)';
end

f_ = K(1); Da_ = K(2); Dep_ = K(3); Dpp_ = K(4); fw_ = K(5);
cosang  = bvecs * dq';                              % [Ndwi x Nq]
Kmat    = f_         * exp(-bvals' .* (Da_*cosang.^2)) + ...
          (1-f_-fw_) * exp(-bvals' .* (Dpp_ + (Dep_-Dpp_)*cosang.^2)) + ...
          fw_        * exp(-bvals' * C.D_FW);
S_exact = zeros(NCOND, Ndwi);
for ic = 1:NCOND
    w = fodf_gt(:,ic) / sum(fodf_gt(:,ic));         % weights, sum to 1
    S_exact(ic,:) = (Kmat * w)';
end

fprintf('Step 4: noise-free signal\n');
fprintf('   shell mean over the 90 acquired directions of each shell:\n');
fprintf('        b   ');  fprintf('%10s', 'single', '15 deg', '45 deg', '60 deg'); fprintf('\n');
for s = shells
    fprintf('     %4g   ', s); fprintf('%10.4f', mean(S_clean(:, bvals==s), 2)); fprintf('\n');
end
fprintf('   These differ by about 1%% across conditions, and that is worth explaining,\n');
fprintf('   because SMI estimates the kernel from rotational invariants and the l = 0\n');
fprintf('   invariant is supposed to be orientation-blind. It is -- but only as an\n');
fprintf('   integral over the whole sphere. Ninety directions do not integrate exactly,\n');
fprintf('   and the residue depends on where the fibres happen to sit. Redo the same\n');
fprintf('   average over the dense quadrature grid and the dependence disappears:\n');

Nq   = size(dq,1);
sm_q = zeros(NCOND, numel(shells));
for i = 1:numel(shells)
    bq = shells(i)*ones(1,Nq);
    for ic = 1:NCOND
        sq = H.signal(plm_gt(ic,:), [K 1 1], bq, ones(1,Nq), zeros(1,Nq), ...
                      dq, C.LMAX_GT, C.CS_PHASE, C.D_FW);
        sm_q(ic,i) = mean(sq);
    end
end
fprintf('        b   ');  fprintf('%10s', 'single', '15 deg', '45 deg', '60 deg'); fprintf('\n');
for i = 1:numel(shells)
    fprintf('     %4g   ', shells(i)); fprintf('%10.6f', sm_q(:,i)); fprintf('\n');
end

% CHECK 0. The dense spherical mean must not depend on the fibre configuration
% at all. If it did, the kernel SMI recovers would depend on the fODF, and the
% whole rotational-invariant approach would be unsound.
% The tolerance is quadrature-limited, not exact: NDIR_Q directions integrate
% the sphere to about 1e-6 here, so a tighter bound would only be testing the
% grid. Raise NDIR_Q in mc_config.m and this residual falls with it.
e_sm = max(max(sm_q, [], 1) - min(sm_q, [], 1));
fprintf('   CHECK spherical mean is orientation-blind  max spread = %.2e   %s\n', ...
        e_sm, VERDICT{1+(e_sm < 1e-5)});

% CHECK 1. S0 must be exactly 1: sigma is set to 1/SNR below, and that is only
% the requested SNR if S0 = 1.
e_s0 = max(abs(S_clean(:,1) - 1));
fprintf('   CHECK S(b=0) == 1                    max|err| = %.2e   %s\n', ...
        e_s0, VERDICT{1+(e_s0 < 1e-12)});

% CHECK 2. Harmonics against direct convolution. These do NOT agree to machine
% precision, and should not: the harmonic route is band limited at LMAX_GT and
% the direct sum is not. The residual below IS the band-limiting error of the
% ground truth, and Step 5 puts it next to the noise.
e_sh = max(abs(S_clean(:) - S_exact(:)));
fprintf('   CHECK harmonics vs direct convolution max|err| = %.2e\n', e_sh);
fprintf('         (this is band limiting, not an error -- it falls with LMAX_GT)\n\n');

%% Step 5 -- Rician noise
% Complex Gaussian noise is added to a real signal and the magnitude taken,
% which is exactly Rician:
%
%  S_noisy = sqrt( (S + sigma*n1)^2 + (sigma*n2)^2 ),   n1, n2 ~ N(0,1)
%
% S0 = 1, so |sigma = 1/SNR|. The seed comes from |mc_config.m|, so this
% walkthrough and the pipeline draw the same noise.

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
fprintf('   the truth''s band-limiting error is %.0fx smaller than sigma, so it is\n', ...
        sigma/e_sh);
fprintf('   not what limits any result here\n');

% CHECK 1. Recover sigma from the simulated data rather than trusting the value
% typed in. This is an independent read of the noise level.
resid = S_noisy - S_rep;
s_hat = std(resid(:));
rel   = abs(s_hat - sigma)/sigma;
fprintf('   CHECK recovered sigma       %.5f vs %.5f, %.1f%% off   %s\n', ...
        s_hat, sigma, 100*rel, VERDICT{1+(rel < 0.10)});

% CHECK 2. A magnitude is strictly positive, which is why SMI is told
% NoiseBias = 'Rician' in Step 6. A negative value would mean the noise model
% is not the one being claimed.
mn = min(S_noisy(:));
fprintf('   CHECK signal strictly positive   min = %.4f            %s\n', ...
        mn, VERDICT{1+(mn > 0)});
hi = (bvals == max(bvals));
fprintf('   Rician floor at b = %g: mean %.4f against a noise-free %.4f (upward bias)\n\n', ...
        max(bvals), mean(mean(S_noisy(:,hi))), mean(mean(S_rep(:,hi))));

%% Step 6 -- SMI.fit
% The real toolbox at its shipped defaults, with the constrained deconvolution
% turned on. Two options are not defaults and both matter:
%
% * |CS_phase = 0|. At SMI's default of 1 the spherical harmonic basis differs
%   from MRtrix's by |(-1)^m|, which for even l is a 180 degree rotation about
%   z of every fODF -- 71.5 degrees of peak error, verified against MRtrix's
%   own |sh2peaks| in |check_mrtrix_basis.sh|. At 0 the two bases are identical.
% * |fODF_regularization.flag_nonneg = 1|. Off by default in the toolbox; this
%   is the arm being studied. |lambda_nonneg| is left at its shipped value and
%   printed back below as confirmation.

options = struct();
options.b     = bvals;
options.dirs  = bvecs;
options.sigma = sigma*ones(GRID);
options.mask  = true(GRID);
options.compartments  = {'IAS','EAS','FW'};
options.NoiseBias     = 'Rician';
options.Lmax          = [0 C.LMAX_FIT C.LMAX_FIT C.LMAX_FIT];
options.CS_phase      = C.CS_PHASE;
options.D_FW          = C.D_FW;
options.flag_fit_fODF = 1;
options.fODF_regularization = struct('flag_nonneg', 1, 'lambda_tikhonov', 0.3);

fprintf('Step 6: SMI.fit on %d voxels ...\n', NVOX);
t0  = tic;
out = SMI.fit(dwi, options);
fprintf('   done in %.1f s   (lambda_nonneg = %g, lambda_tikhonov = %g)\n', ...
        toc(t0), out.fODF_regularization.lambda_nonneg, ...
        out.fODF_regularization.lambda_tikhonov);

kern = reshape(out.kernel, [NVOX size(out.kernel,4)]);
knm  = {'f','Da','Depar','Deperp','fw'};
fprintf('   kernel recovery, median over all %d voxels:\n', NVOX);
fprintf('        param      true    median      bias\n');
for j = 1:5
    col = kern(:,j); col = col(isfinite(col));
    md  = median(col);
    fprintf('     %9s  %7.3f   %7.3f   %+7.3f\n', knm{j}, K(j), md, md - K(j));
end
fprintf(['   Read that table before reading the fODF results. The kernel is NOT recovered\n' ...
         '   exactly -- a 3-shell linear-encoding protocol constrains it only loosely, which\n' ...
         '   is what README.md "Useful tips" warns about, and SMI''s estimator is a regression\n' ...
         '   trained on a prior, so biased parameters get pulled toward that prior.\n' ...
         '   The fODF below is then deconvolved with this ESTIMATED kernel, not the true one.\n' ...
         '   That is deliberate: it is what happens on real data, where nobody has the truth.\n' ...
         '   Deconvolving with the true kernel instead is a one-line change and a fair\n' ...
         '   experiment to run -- it separates kernel error from deconvolution error.\n']);

% Convert plm to spherical harmonic coefficients in the same convention the
% ground truth was written in, so the two are directly comparable.
L6  = repelem(0:2:C.LMAX_FIT, 2*(0:2:C.LMAX_FIT)+1)';
sc6 = sqrt((2*L6+1)/(4*pi))';
plm = reshape(out.plm, [NVOX numel(L6)-1]);
sh  = [ones(NVOX,1) plm] .* repmat(sc6, NVOX, 1);
nonfinite = sum(~isfinite(sh(:)));
sh(~isfinite(sh)) = 0;

% CHECK 1. The constrained deconvolution must converge in every voxel. A voxel
% that hit the iteration cap has not satisfied the constraint, so its fODF is
% not the thing the method promises.
conv = out.fODF_regularization.flag_converged(:);
fprintf('   CHECK deconvolution converged   %d / %d voxels         %s\n', ...
        sum(conv == 1), NVOX, VERDICT{1+all(conv == 1)});

% CHECK 2. No NaN reached the fODF: a NaN in an SH volume breaks downstream
% tractography silently.
fprintf('   CHECK fODF all finite           %d non-finite          %s\n', ...
        nonfinite, VERDICT{1+(nonfinite == 0)});

% CHECK 3. The l = 0 coefficient must still be the fixed isotropic term, i.e.
% the p_00 = 1 convention survived the fit.
e_l0 = max(abs(sh(:,1) - 1/sqrt(4*pi)));
fprintf('   CHECK p_00 == 1 convention held max|err| = %.2e   %s\n\n', ...
        e_l0, VERDICT{1+(e_l0 < 1e-12)});

%% Step 7 -- peaks and angular error
% Did it work? The only questions a tractography algorithm asks of an fODF are
% "how many fibres, and pointing where", so those are what get scored.
%
% Peaks are found by evaluating the fODF on a dense direction set and keeping
% every direction not smaller than any neighbour within |PEAK_NBR| degrees,
% then keeping those whose *anisotropic* amplitude -- amplitude minus the
% 1/(4*pi) floor -- is at least |PEAK_REL| of the largest. Subtracting the
% floor matters: without it a nearly isotropic fODF looks like it has many
% strong peaks.
%
% *This is the walkthrough's own peak finder, not the report's.* The report
% puts every arm through MRtrix's |sh2peaks| so peak extraction is identical
% across methods. Step 8 checks the two against each other.

PEAK_NBR = 12;      % degrees, neighbourhood for the local-maximum test
PEAK_REL = 0.30;    % keep peaks at >= 30% of the largest anisotropic amplitude
de  = H.dirs(1500);
Ye  = SMI.get_even_SH(de, C.LMAX_FIT, C.CS_PHASE);
nbr = (de*de') > cosd(PEAK_NBR);
nbr(1:size(de,1)+1:end) = true;               % include self, so >= means "is the max"

A_fit = sh    * Ye' - 1/(4*pi);               % anisotropic amplitude, per voxel
A_gt  = sh_gt * Ye' - 1/(4*pi);

fprintf('Step 7: peaks (%d directions, %d deg neighbourhood, %.0f%% threshold)\n', ...
        size(de,1), PEAK_NBR, 100*PEAK_REL);
fprintf('   condition   correct count   median angular error   spurious peaks\n');

res_frac = zeros(1,NCOND); ang_med = zeros(1,NCOND); spur = zeros(1,NCOND);
for ic = 1:NCOND
    rows   = find(cond_id == ic);
    ntrue  = numel(axes_gt{ic});
    nfound = zeros(numel(rows),1);
    aerr   = nan(numel(rows),1);
    for r = 1:numel(rows)
        a  = A_fit(rows(r),:)';
        lm = find(a > 0 & a >= max(a .* nbr, [], 1)');
        if isempty(lm), continue; end
        lm = lm(a(lm) >= PEAK_REL*max(a(lm)));
        [~, o] = sort(a(lm), 'descend');
        lm = lm(o);
        P  = de(lm,:);
        sel = true(size(P,1),1);              % merge antipodal duplicates
        for i = 1:size(P,1)
            if ~sel(i), continue; end
            dup = abs(P*P(i,:)') > cosd(PEAK_NBR); dup(i) = false;
            sel(dup) = false;
        end
        P = P(sel,:);
        nfound(r) = size(P,1);
        d = zeros(1,ntrue);
        for k = 1:ntrue, d(k) = acosd(min(abs(P(1,:)*axes_gt{ic}{k}'), 1)); end
        aerr(r) = min(d);
    end
    res_frac(ic) = 100*mean(nfound == ntrue);
    ang_med(ic)  = median(aerr(isfinite(aerr)));
    spur(ic)     = mean(max(nfound - ntrue, 0));
    fprintf('   %6d deg   %11.1f%%   %16.2f deg   %14.3f\n', ...
            C.ANGLES(ic), res_frac(ic), ang_med(ic), spur(ic));
end

% The ceiling. Run the identical finder on the band-limited ground truth, with
% no noise and no fitting. Nothing above can beat this, and two things are read
% off it.
%
% * *The 15 degree crossing is not resolvable, by construction.* The truth
%   itself comes back as a single peak. The 0.0% above is therefore the correct
%   answer, not a failure -- two Watson populations with kappa = 16 only 15
%   degrees apart merge into one lobe long before any deconvolution is involved.
% * *The finder's own resolution floor.* A single fibre whose truth is "1.7 deg
%   off" is measuring the spacing of the evaluation grid, not an error. The
%   report avoids this floor by using |sh2peaks|, which refines each peak with
%   a Newton search on the continuous expansion instead of snapping to a grid.

fprintf('   the ceiling -- the same finder on the band-limited truth, no noise, no fit:\n');
ceil_n = zeros(1,NCOND);
for ic = 1:NCOND
    a  = A_gt(ic,:)';
    lm = find(a > 0 & a >= max(a .* nbr, [], 1)');
    lm = lm(a(lm) >= PEAK_REL*max(a(lm)));
    [~, o] = sort(a(lm), 'descend'); lm = lm(o);
    P = de(lm,:);
    sel = true(size(P,1),1);
    for i = 1:size(P,1)
        if ~sel(i), continue; end
        dup = abs(P*P(i,:)') > cosd(PEAK_NBR); dup(i) = false; sel(dup) = false;
    end
    P = P(sel,:);
    ntrue     = numel(axes_gt{ic});
    ceil_n(ic) = size(P,1);
    d = zeros(1,ntrue);
    for k = 1:ntrue, d(k) = acosd(min(abs(P(1,:)*axes_gt{ic}{k}'), 1)); end
    if ceil_n(ic) == ntrue, note = 'resolvable';
    else,                   note = 'NOT resolvable even from the truth';
    end
    fprintf('   %6d deg   %d peak(s) of a true %d, largest %.2f deg off   %s\n', ...
            C.ANGLES(ic), ceil_n(ic), ntrue, min(d), note);
end

% CHECK. Where the truth itself cannot be resolved, the fit must not claim to
% resolve it. A method that "finds" two fibres where the noise-free truth has
% one lobe is inventing them, and that is the failure mode worth catching.
unres = false; ok_unres = true;
for ic = 1:NCOND
    if ceil_n(ic) < numel(axes_gt{ic})
        unres = true;
        ok_unres = ok_unres && (res_frac(ic) < 5);
    end
end
if unres
    fprintf('   CHECK fit does not resolve what the truth cannot        %s\n', ...
            VERDICT{1+ok_unres});
end

% CHECK. The constrained deconvolution's whole claim is that it does not
% manufacture peaks. At SNR 50 the spurious count should be flat zero; the
% report finds it stays at 0.000 down to SNR 10 and only rises at SNR 5.
fprintf('   CHECK no spurious peaks            max = %.3f            %s\n', ...
        max(spur), VERDICT{1+(max(spur) < 0.01)});
fprintf('\n');

%% Step 8 -- cross-check the peak finder against MRtrix
% The report does not use the finder above; it runs every arm through MRtrix's
% |sh2peaks| so that peak extraction is identical across methods and is not
% something this package implements. If MRtrix3 is installed the two are
% compared here on the same fODFs. If it is not, this section is skipped and
% nothing else depends on it.

if system('command -v sh2peaks > /dev/null 2>&1') ~= 0
    fprintf('Step 8: sh2peaks not found, skipping the cross-check.\n');
    fprintf('        Install MRtrix3 (apt-get install -y mrtrix3) and rerun.\n\n');
else
    MR   = mrtrix_io();
    mdir = fullfile(pkgdir, 'mrtrix');
    if ~exist(mdir, 'dir'), mkdir(mdir); end
    fod_f = fullfile(mdir, 'walkthrough_smifod');
    pk_f  = fullfile(mdir, 'walkthrough_peaks');
    MR.write(fod_f, reshape(sh, [GRID numel(L6)]));
    st = system(sprintf('sh2peaks -quiet -force -num 4 "%s.mih" "%s.mih"', fod_f, pk_f));
    if st ~= 0
        fprintf('Step 8: sh2peaks exited %d, skipping.\n\n', st);
    else
        pk = reshape(MR.read([pk_f '.mih']), [NVOX 12]);
        fprintf('Step 8: correct fibre count, sh2peaks vs this walkthrough\n');
        fprintf('   condition     sh2peaks   walkthrough\n');
        for ic = 1:NCOND
            rows  = find(cond_id == ic);
            ntrue = numel(axes_gt{ic});
            nmr   = zeros(numel(rows),1);
            for r = 1:numel(rows)
                v   = reshape(pk(rows(r),:), 3, 4)';
                nrm = sqrt(sum(v.^2, 2));
                nrm = nrm(isfinite(nrm) & nrm > 0);
                if isempty(nrm), continue; end
                nmr(r) = sum(nrm >= PEAK_REL*max(nrm));
            end
            fprintf('   %6d deg   %9.1f%%   %10.1f%%\n', ...
                    C.ANGLES(ic), 100*mean(nmr == ntrue), res_frac(ic));
        end
        fprintf(['   These need not match exactly. sh2peaks runs a Newton search on the\n' ...
                 '   continuous SH expansion and thresholds raw amplitude; the finder above\n' ...
                 '   searches a fixed grid and thresholds anisotropic amplitude. They should\n' ...
                 '   agree on the trend. Every number in the report comes from sh2peaks.\n\n']);
    end
end

%% Step 9 -- from here to the published run
% Everything above is one SNR at |NREP| realisations. The published campaign is
% the same code at |NREP = 10000| across six noise levels, driven by
% |gen_montecarlo.m| rather than by hand:
%
%  cd deconv_comparison
%  octave --eval "run('oct_path.m'); gen_montecarlo(50, 10000, 'snr50')"
%  ./run_mrtrix.sh fit snr50            % CSD, MSMT-CSD and every peak, in MRtrix
%  python3 tables.py nf snr50:50 ...    % the tables in Reports/deconv_tables.md
%
% or |./run_all.sh| for all of it (|NREP=200 ./run_all.sh| for a quick pass).
%
% *What this walkthrough leaves out on purpose.* The comparison arms. CSD and
% MSMT-CSD are MRtrix3 binaries driven by |run_mrtrix.sh|, and the scoring that
% puts all three side by side is |score_mrtrix.py|. This file is the SMI arm
% alone, taken apart.
%
% *What to change first if you want to probe it.* |NREP| and |SNR| at the top
% of this file. Then |mc_config.m|: |KAPPA| sets fibre dispersion, |K_WM| the
% tissue, |ANGLES| the crossings, |LMAX_FIT| what the methods may represent.
% Changing them there changes them for the pipeline too, which is the point of
% that file existing.

fprintf('=== walkthrough complete ===\n');
fprintf('SNR %g, %d reps per condition, correct fibre count:\n', SNR, NREP);
for ic = 1:NCOND
    fprintf('   %2d deg  %5.1f%%\n', C.ANGLES(ic), res_frac(ic));
end
fprintf('\n');
