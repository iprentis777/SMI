function C = mc_config()
% C = mc_config()
%
% ONE definition of the Monte Carlo experiment. gen_montecarlo.m and
% sweep_nonneg.m both read their constants from here, so the main run and the
% regularizer sweep cannot disagree about what is being simulated. Anything
% added later that narrates or re-runs the experiment should read them from
% here too rather than restating them.
%
% Constants
%   LMAX_FIT   angular order every method deconvolves at
%   LMAX_GT    angular order of the ground truth (8 is the SM kernel's
%              ceiling, SMI.m:2470-2480), deliberately ABOVE LMAX_FIT so no
%              method is handed a truth it can represent exactly
%   CS_PHASE   0 == MRtrix's SH basis exactly. At SMI's default of 1 the two
%              bases differ by (-1)^m, a 180 degree rotation about z of every
%              fODF; check_mrtrix_basis.sh measures that against MRtrix itself
%   D_FW       free water diffusivity, um^2/ms
%   KAPPA      Watson concentration of each fibre population. Finite, not a
%              delta: a response estimated from real white matter absorbs
%              fibre dispersion, and a delta truth would hand every method a
%              mismatch that does not exist in practice
%   ANGLES     crossing angles in degrees; 0 means a single fibre. NOTE that at
%              KAPPA = 16 a 30 degree crossing does not separate even in the
%              noise-free truth -- measured at Lmax 8 it needs KAPPA >= 48.
%              45 and 60 separate at every KAPPA from 8 upward
%   K_WM       [f Da Depar Deperp fw] of the white matter kernel
%   AXIS1      the first fibre axis, fixed and off every coordinate plane
%   NDIR_Q     quadrature directions used to project a sampled fODF onto plm
%   SEED_MC, SEED_SWEEP   RNG seeds, one per script, so the two are independent
%
% Helpers
%   C.condition_axes(ic)  the true fibre axes of condition ic, as a cell array
%   C.pick_grid(N)        a 3D grid holding exactly N voxels, no singleton dim
%   C.rotate_about(n,deg) a unit vector `deg` degrees from n, in a fixed plane

C = struct();
C.PROTOCOL   = 'hcp_real_3shell.txt';
C.LMAX_FIT   = 6;
C.LMAX_GT    = 8;
C.CS_PHASE   = 0;
C.D_FW       = 3;
C.KAPPA      = 16;
C.ANGLES     = [0 30 45 60];
C.K_WM       = [0.60 2.0 2.0 0.50 0.02];
C.NDIR_Q     = 3000;
C.SEED_MC    = 31415;
C.SEED_SWEEP = 2718;

n1 = [0.30 -0.50 0.81];
C.AXIS1 = n1/norm(n1);

C.pick_grid     = @pick_grid;
C.rotate_about  = @rotate_about;
C.condition_axes = @(ic) condition_axes(C.AXIS1, C.ANGLES, ic);
C.load_protocol  = @() load_protocol(C.PROTOCOL);
end

% =====================================================================
function [bvals, bvecs] = load_protocol(fname)
% [bvals, bvecs] = load_protocol(fname)
%
% The acquisition, read from protocol/<fname>. b comes back in ms/um^2 as a
% row; bvecs is [N x 3] and unit.
%
% Every arm reads the protocol through here so there is one definition of what
% was acquired, and so the normalisation warning below cannot be bypassed by
% loading the file directly.
here = fileparts(mfilename('fullpath'));
txt  = fileread(fullfile(here, 'protocol', fname));
cols = textscan(txt, '%f %f %f %f', 'CommentStyle', '%');
bvals = cols{1}(:)';
bvecs = [cols{2} cols{3} cols{4}];

% Real .bvec files are written at finite precision, so the directions are not
% quite unit. This is NOT harmless, and it is not silently repaired: anything
% that treats g(3) as cos(theta) without normalising inherits the error, and at
% Lmax 8 an error of 1e-6 in |g| is enough to break the zonal-response identity
% by 5e-7 -- six orders of magnitude worse than machine precision. Warn with the
% measured number, then correct it.
nrm   = sqrt(sum(bvecs.^2, 2));
e_raw = max(abs(nrm - 1));
if e_raw > 1e-12
    fprintf(2, ['\n' repmat('!',1,72) '\n']);
    fprintf(2, 'WARNING  %s: gradient directions are not unit vectors.\n', fname);
    fprintf(2, '         max | |g| - 1 | = %.3e over %d volumes.\n', e_raw, numel(bvals));
    fprintf(2, '         They are being normalised. Left uncorrected this breaks any\n');
    fprintf(2, '         calculation that reads g(3) as cos(theta) -- measured at Lmax 8,\n');
    fprintf(2, '         the zonal response vs forward model identity degrades from\n');
    fprintf(2, '         ~1e-15 to ~5e-7. Check the source of the .bvec if this is large.\n');
    fprintf(2, [repmat('!',1,72) '\n\n']);
end
bvecs = bvecs ./ repmat(nrm, 1, 3);
end

% =====================================================================
function axes_ = condition_axes(n1, ANGLES, ic)
% The true fibre axes of condition ic: one axis for the single fibre case,
% two equal populations otherwise.
if ANGLES(ic) == 0
    axes_ = {n1};
else
    axes_ = {n1, rotate_about(n1, ANGLES(ic))};
end
end

% =====================================================================
function m = rotate_about(n, deg)
% A unit vector at `deg` degrees from n, in an arbitrary but FIXED plane, so
% every script that asks for the same angle gets the same pair of axes.
n = n(:)'/norm(n);
t = [0 0 1]; if abs(n*t') > 0.9, t = [1 0 0]; end
e = t - (t*n')*n; e = e/norm(e);
m = cosd(deg)*n + sind(deg)*e;
m = m/norm(m);
end

% =====================================================================
function G = pick_grid(N)
% A 3D grid holding exactly N voxels with every dimension > 1: SMI.vectorize
% takes a different branch if any spatial dimension is a singleton
% ("README for Claude", section 4). The most cube-like factorisation wins,
% purely so the printed grid looks sensible; nothing depends on the shape.
d = divisors_of(N);
d = d(d > 1 & d < N);
best = [];
for a = d
    m = N/a;
    e = divisors_of(m);
    e = e(e > 1 & e < m);
    for b = e
        c = m/b;
        if c > 1
            cand = sort([a b c]);
            if isempty(best) || (max(cand)-min(cand)) < (max(best)-min(best))
                best = cand;
            end
        end
    end
end
if isempty(best)
    error('mc_config:pick_grid', ...
          'cannot factor %d into three factors > 1 -- pick another NREP', N);
end
G = best;
end

% =====================================================================
function d = divisors_of(n)
d = 1:floor(sqrt(n));
d = d(mod(n,d) == 0);
d = unique([d n./d]);
end
