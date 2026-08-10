function nfail = test_csd_arms(neg_sweep)
% nfail = test_csd_arms([neg_sweep])
%
% Exercise the CSD arms on their own, in about half a minute, without waiting on
% SMI.fit.
%
%   octave-cli --no-gui -q test_csd_arms.m
%   >> test_csd_arms                 % the standard checks
%   >> test_csd_arms(true)           % plus the -neg_lambda sensitivity sweep
%
% smi_manuscript_60deg.m runs all three arms, but its SMI arm costs 42 fits and
% hours at the manuscript settings, so it is the wrong instrument for "does the
% MRtrix side still work". This file builds the same noise-free signal from the
% same kernel and protocol, hands it to dwi2fod, and scores the peaks -- the
% whole MRtrix path, none of the cost.
%
% WHAT IT ASSERTS, AND WHY EACH ONE IS HERE
%
% 1. *MRtrix reads what we wrote.* The .mif writer is ours; mrinfo is the only
%    thing that can confirm it.
% 2. *The shells survive the round trip*, and MRtrix's own per-shell b values
%    are what the response is evaluated at. This protocol's b values jitter.
% 3. *SSST-CSD reproduces the band-limited truth* on a noise-free crossing. This
%    is the reference: if it fails, the basis, the response convention or the
%    peak finder is broken, and no other number here means anything.
% 4. *MSMT-CSD separates the crossing as well as SSST-CSD does.* This is the
%    regression test for the bug that made it into the repository once already:
%    dwi2fod csd and dwi2fod msmt_csd do NOT ship comparable regularisation, and
%    running both "at their defaults" compares a constrained arm against an
%    effectively unconstrained one. At MRtrix's -neg_lambda default of 1e-10 the
%    MSMT fODF comes back blunt and under-separates a 60 degree crossing by
%    ~12 degrees. See the Configuration block of smi_manuscript_60deg.m.
% 5. *90 degrees works.* A method that cannot separate an orthogonal crossing
%    noise-free is misconfigured, whatever it does at 60. This is the check that
%    would have caught the bug immediately.
%
% Pass `true` to also sweep -neg_lambda and -norm_lambda and print the
% separation for each. That sweep is worth reading before quoting any MSMT
% number: the answer moves by 16 degrees across the range, and the norm_lambda
% value in the manuscript file is the least justified constant in this package.

if nargin < 1 || isempty(neg_sweep), neg_sweep = false; end
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'oct_path.m'));

MC = mc_config(); H = fODF_modulation_helpers();
RH = SMI_response_helpers(); MR = mrtrix_io(); PK = fODF_peak_score();
VERDICT = {'** FAILED **', 'ok'};
nfail = 0;

% ---- the same configuration smi_manuscript_60deg.m uses
K = [0.60 2.0 2.0 0.50 0.02]; KAPPA = 16; D_FW = 3; D_GM = 0.8;
LGT = 8; CS = 0; Lf = 6; B0_SNAP = 0.05;
NEG_LAMBDA = 1; NORM_LAMBDA = 1e-3;      % must match the manuscript file
AXIS1 = [0.30 -0.50 0.81]; AXIS1 = AXIS1/norm(AXIS1);

D = fullfile(here, 'mrtrix');
if ~exist(D, 'dir'), mkdir(D); end

fprintf('\n=== CSD arms, standalone (no SMI.fit) ===\n');
[bvals, bvecs] = MC.load_protocol_file('hcp_real_3shell.txt');
bvals(bvals < B0_SNAP) = 0;
dq = H.dirs(3000); de = H.dirs(1500);
ctx = PK.setup(de, 12, 0.30);

% Watson zonal coefficients, for the dispersion-matched response
plm_w = H.mixture_plm(H.watson(dq, [0 0 1], KAPPA), dq, LGT, CS);
Mv = []; for il = 0:2:LGT, Mv = [Mv, -il:il]; end %#ok<AGROW>
sh_w = [1; plm_w(:)]; PL = sh_w(Mv(:) == 0);

G = [3 3 3]; NV = prod(G);
f_mask = fullfile(D, 'tca_mask');
MR.write([f_mask '.mif'], ones(G), struct('datatype', 'UInt8'));

angles = [60 75 90];
sep_ssst = nan(size(angles)); sep_msmt = nan(size(angles)); sep_true = nan(size(angles));

for ia = 1:numel(angles)
    ang = angles(ia);
    ax2 = MC.rotate_about(AXIS1, ang);
    fodf = H.watson(dq, AXIS1, KAPPA) + H.watson(dq, ax2, KAPPA);
    plm  = H.mixture_plm(fodf, dq, LGT, CS);
    S    = H.signal(plm, [K 1 1], bvals, ones(1,numel(bvals)), ...
                    zeros(1,numel(bvals)), bvecs, LGT, CS, D_FW);

    f_dwi = fullfile(D, sprintf('tca_%d', ang));
    MR.write([f_dwi '.mif'], reshape(repmat(S(:)', NV, 1), [G numel(bvals)]), ...
             struct('grad', [bvecs bvals(:)*1000]));

    if ia == 1
        [st, txt] = system(sprintf('mrinfo -size "%s.mif" 2>&1', f_dwi));
        sz = sscanf(txt, '%f')';
        ok = (st == 0) && numel(sz) == 4 && all(sz == [G numel(bvals)]);
        nfail = nfail + ~ok;
        fprintf('   CHECK mrinfo reads back %s   %s\n', mat2str(sz), VERDICT{1+ok});
        if ~ok
            fprintf(2, '%s\n', txt);
            error(['MRtrix could not read the image this file wrote. ' ...
                   'Is mrtrix3 on the PATH?']);
        end
        [~, txt] = system(sprintf('mrinfo -shell_sizes "%s.mif" 2>/dev/null', f_dwi));
        ns = sscanf(txt, '%f')';
        ok = isequal(ns, [18 90 90 90]);
        nfail = nfail + ~ok;
        fprintf('   CHECK shells come back %s   %s\n', mat2str(ns), VERDICT{1+ok});
    end

    [~, txt] = system(sprintf('mrinfo -shell_bvalues "%s.mif" 2>/dev/null', f_dwi));
    b_mr = sscanf(txt, '%f')';

    r_wm  = RH.zh(K, b_mr/1000, Lf, D_FW) .* repmat(PL(1:Lf/2+1)', numel(b_mr), 1);
    r_gm  = exp(-(b_mr(:)/1000)*D_GM)*sqrt(4*pi);
    r_csf = exp(-(b_mr(:)/1000)*D_FW)*sqrt(4*pi);
    Rw = fullfile(D,'tca_wm.txt');  RH.write_response(Rw,  r_wm);
    Rg = fullfile(D,'tca_gm.txt');  RH.write_response(Rg,  r_gm);
    Rc = fullfile(D,'tca_csf.txt'); RH.write_response(Rc,  r_csf);
    Rb = fullfile(D,'tca_b3.txt');  RH.write_response(Rb,  r_wm(end,:));

    f_b3 = fullfile(D, sprintf('tca_%d_b3', ang));
    system(sprintf('dwiextract "%s.mif" -shells 0,%g "%s.mif" -force -quiet 2>&1', ...
                   f_dwi, b_mr(end), f_b3));
    system(sprintf('dwi2fod csd "%s.mif" "%s" "%s/tca_s.mif" -mask "%s.mif" -force -quiet 2>&1', ...
                   f_b3, Rb, D, f_mask));
    system(sprintf(['dwi2fod msmt_csd "%s.mif" "%s" "%s/tca_m.mif" "%s" "%s/tca_mg.mif" ' ...
                    '"%s" "%s/tca_mc.mif" -mask "%s.mif" -neg_lambda %g -norm_lambda %g ' ...
                    '-force -quiet 2>&1'], ...
                   f_dwi, Rw, D, Rg, D, Rc, D, f_mask, NEG_LAMBDA, NORM_LAMBDA));

    L8 = repelem(0:2:LGT, 2*(0:2:LGT)+1)';
    sh_gt = ([1; plm(:)] .* sqrt((2*L8+1)/(4*pi)))';
    nc = (Lf/2+1)*(Lf+1);

    sep_true(ia) = sep_of(sh_gt(1:nc), de, ctx, PK, Lf, CS);
    sep_ssst(ia) = sep_of(read1(MR, sprintf('%s/tca_s.mif', D)), de, ctx, PK, Lf, CS);
    sep_msmt(ia) = sep_of(read1(MR, sprintf('%s/tca_m.mif', D)), de, ctx, PK, Lf, CS);
end

fprintf('\n   peak separation, noise free, Lmax %d, dispersion-matched response\n', Lf);
fprintf('   msmt_csd -neg_lambda %g -norm_lambda %g;  csd constraint strength 1\n', ...
        NEG_LAMBDA, NORM_LAMBDA);
fprintf('     true    band-limited truth    SSST-CSD    MSMT-CSD\n');
for ia = 1:numel(angles)
    fprintf('     %3d deg        %6.2f          %6.2f      %6.2f\n', ...
            angles(ia), sep_true(ia), sep_ssst(ia), sep_msmt(ia));
end

% SSST is the reference, and it is scored against the TRUE angle rather than
% against the band-limited truth column.
%
% The first version of this check compared the two and failed at 75 degrees,
% where the truncated truth separates at 70.99 while both CSD arms give 73.95 --
% i.e. the arms beat the "ceiling". That is not an error and the check was
% wrong: with a DISPERSION-MATCHED response the fODF the arms are asked to
% recover is a pair of deltas, not the truncated Watson mixture, so they are not
% bounded by how well a truncated Watson separates. The band-limited column is
% still printed, because the gap between it and the arms is worth seeing, but it
% is not the bound.
e = max(abs(sep_ssst - angles));
ok = e < 2.6;                      % the direction grid's own resolution floor
nfail = nfail + ~ok;
fprintf('\n   CHECK SSST-CSD recovers the true angle            max|err| = %.2f deg   %s\n', ...
        e, VERDICT{1+ok});

% The regression test. MSMT must not be systematically blunter than SSST.
d = max(sep_ssst - sep_msmt);
ok = d < 2.0;
nfail = nfail + ~ok;
fprintf('   CHECK MSMT-CSD separates as well as SSST-CSD      max gap  = %.2f deg   %s\n', ...
        d, VERDICT{1+ok});
if ~ok
    fprintf(2, ['        MSMT is under-separating. The usual cause is -neg_lambda:\n' ...
                '        MRtrix defaults it to 1e-10, which leaves msmt_csd effectively\n' ...
                '        unconstrained while csd is constrained at strength 1. Run\n' ...
                '        test_csd_arms(true) to see the sensitivity.\n']);
end

% 90 degrees is the sanity check against the literature.
i90 = find(angles == 90, 1);
if ~isempty(i90)
    ok = (sep_msmt(i90) > 85) && (sep_ssst(i90) > 85);
    nfail = nfail + ~ok;
    fprintf('   CHECK both arms separate a 90 deg crossing        %s\n', VERDICT{1+ok});
end

% ---- optional: the sensitivity that makes these knobs worth reporting
if neg_sweep
    fprintf('\n   -neg_lambda / -norm_lambda sensitivity at a true 60 deg crossing\n');
    fprintf('   (band-limited truth separates at %.2f deg)\n', sep_true(1));
    fprintf('     neg_lambda   norm_lambda   separation\n');
    ang = 60; f_dwi = fullfile(D, 'tca_60');
    [~, txt] = system(sprintf('mrinfo -shell_bvalues "%s.mif" 2>/dev/null', f_dwi));
    b_mr = sscanf(txt, '%f')';
    r_wm = RH.zh(K, b_mr/1000, Lf, D_FW) .* repmat(PL(1:Lf/2+1)', numel(b_mr), 1);
    RH.write_response(fullfile(D,'tca_wm.txt'), r_wm);
    for nl = [1e-10 1e-3 1]
        for gl = [1e-10 1e-3 1]
            system(sprintf(['dwi2fod msmt_csd "%s.mif" "%s/tca_wm.txt" "%s/tca_q.mif" ' ...
                            '"%s/tca_gm.txt" "%s/tca_qg.mif" "%s/tca_csf.txt" "%s/tca_qc.mif" ' ...
                            '-mask "%s.mif" -neg_lambda %g -norm_lambda %g -force -quiet 2>&1'], ...
                           f_dwi, D, D, D, D, D, D, f_mask, gl, nl));
            fprintf('     %10g   %11g   %8.2f\n', gl, nl, ...
                    sep_of(read1(MR, sprintf('%s/tca_q.mif', D)), de, ctx, PK, Lf, CS));
        end
    end
    fprintf(['     MRtrix ships neg_lambda = norm_lambda = 1e-10. Report whichever\n' ...
             '     values you use: the answer moves by ~16 deg across this table.\n']);
end

fprintf('\n');
if nfail == 0
    fprintf('All CSD checks passed.\n\n');
else
    fprintf('%d CSD check(s) FAILED.\n\n', nfail);
end
end

% =====================================================================
function v = read1(MR, fname)
% The SH coefficients of the first voxel of an MRtrix image.
V = MR.read(fname);
v = reshape(V, [], size(V,4));
v = v(1,:);
end

% =====================================================================
function s = sep_of(v, de, ctx, PK, Lf, CS)
% Angle between the two strongest peaks of an fODF given as SH coefficients.
% The isotropic part removed is the voxel's OWN l = 0 term, which is what makes
% this work for an unnormalised MRtrix FOD as well as for an SMI one.
Ye = SMI.get_even_SH(de, Lf, CS);
A  = v(1:size(Ye,2))*Ye' - v(1)/sqrt(4*pi);
P  = PK.peaks(A, ctx);
if size(P,1) >= 2
    s = acosd(min(abs(P(1,:)*P(2,:)'), 1));
else
    s = NaN;
end
end
