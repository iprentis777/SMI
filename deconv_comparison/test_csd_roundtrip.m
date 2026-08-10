function test_csd_roundtrip()
% test_csd_roundtrip -- re-score the Python CSD arms with the Octave peak finder.
%
%   octave-cli --no-gui -q test_csd_roundtrip.m
%   >> test_csd_roundtrip            % from MATLAB or Octave
%
% The notebook exports both the SH coefficients AND the numbers it scored them
% to. That redundancy is the point: this file loads the coefficients, runs the
% OCTAVE peak finder over them, and compares against what Python printed. Two
% independent implementations of the same algorithm, on the same input.
%
% What each outcome means:
%
%   * both agree            -- the scoring is implementation independent, and
%                              the .m file can trust the exported score arrays
%                              without re-deriving them.
%   * counts agree, angles differ slightly
%                           -- the peak grid resolves the same lobes but breaks
%                              a near-tie differently. Expected on a SYMMETRIC
%                              crossing, where the two lobes have near-equal
%                              amplitude; see helpers/fODF_peak_score.m.
%   * counts differ         -- a real disagreement. Look at the SH basis and the
%                              isotropic subtraction first.
%
% It also checks what a silent format change would break: the array shapes, the
% voxel layout, and that the exported noisy signal really is data at the SNRs
% claimed for it.
%
% This is a FUNCTION file, not a script, because Octave cannot call functions
% defined at the end of a script while MATLAB requires them there. Subfunctions
% of a function file work in both -- the same reason helpers/ is laid out the
% way it is ("README for Claude", section 4).

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'oct_path.m'));

PK = fODF_peak_score();
H  = fODF_modulation_helpers();
nfail = 0;

fprintf('\n=== the Python CSD export, re-scored in Octave ===\n');
C = read_csd_export();

NCOND = numel(C.angles);
NSNR  = numel(C.snr_list);
NARM  = numel(C.arms);
NL    = numel(C.lmax_list);

% ---- the peak finder, set up exactly as the notebook set up its own
NDIR_E = 1500;
de  = H.dirs(NDIR_E);
ctx = PK.setup(de, 12, 0.30);
fprintf('\nOctave peak finder: %d directions, 12 deg neighbourhood, 30%% threshold\n', NDIR_E);

% ---- shapes and layout, before any numbers are compared
nfail = nfail + ~check('score array has the documented shape', ...
    isequal(size(C.scores), [NARM NL NSNR NCOND 5]), mat2str(size(C.scores)));
nfail = nfail + ~check('every voxel has a condition and an SNR', ...
    numel(C.cond_id) == C.nvox && numel(C.snr_id) == C.nvox, '');
nfail = nfail + ~check('the grid holds exactly the voxels', ...
    prod(C.grid_all) == C.nvox, sprintf('%s = %d', mat2str(C.grid_all), C.nvox));
nfail = nfail + ~check('no dimension of the grid is a singleton', ...
    all(C.grid_all > 1), 'SMI.vectorize takes a different branch if one is');

% ---- the exported signal really is noisy data at the SNRs claimed
fprintf('\nthe exported signal, noise level recovered per block:\n');
for is = 1:NSNR
    rows = find(C.snr_id == is);
    blk  = C.signal(rows, :);
    if isinf(C.snr_list(is))
        % Every realisation of the noise-free block is the same signal vector.
        d = max(max(abs(blk - repmat(blk(1,:), size(blk,1), 1))));
        nfail = nfail + ~check('SNR inf  block is identical row to row', ...
            d == 0, sprintf('max|row-row1| = %.2e', d));
    else
        % Only b = 0 is a known constant (S0 = 1), so the noise is read there.
        b0    = C.bvals == 0;
        s_hat = std(reshape(blk(:,b0) - 1, [], 1));
        sg    = 1/C.snr_list(is);
        rel   = abs(s_hat - sg)/sg;
        nfail = nfail + ~check(sprintf('SNR %-4g recovered sigma %.5f vs %.5f', ...
            C.snr_list(is), s_hat, sg), rel < 0.15, sprintf('%.1f%% off', 100*rel));
    end
end

% ---- the re-scoring itself
fprintf('\nre-scoring %d arms x %d Lmax x %d SNR:\n', NARM, NL, NSNR);
fprintf('     arm         Lmax   SNR    python err   octave err    python n   octave n\n');
max_derr = 0; max_dspur = 0; max_dres = 0;
for ia = 1:NARM
    for iL = 1:NL
        Lf = C.lmax_list(iL);
        Ye = SMI.get_even_SH(de, Lf, C.cs_phase);
        nc = size(Ye,2);
        for ic = 1:NCOND
            ntrue = numel(C.axes{ic});
            for is = 1:NSNR
                rows = find(C.cond_id == ic & C.snr_id == is);
                % Row 1 is the band-limited truth, as the notebook does it, so
                % the ceiling goes through the same lines as the realisations.
                blk = [C.truth_sh(ic,1:nc); C.sh{ia}{iL}(rows,:)];
                [nf, ae] = PK.score(blk, Ye, C.axes{ic}, ctx);
                nf = nf(2:end); ae = ae(2:end);
                fin = isfinite(ae);

                oct_err = NaN;
                if any(fin), oct_err = mean(ae(fin)); end
                oct_res  = 100*mean(nf == ntrue);
                oct_spur = mean(max(nf - ntrue, 0));

                py_res  = C.scores(ia,iL,is,ic,1);
                py_err  = C.scores(ia,iL,is,ic,2);
                py_spur = C.scores(ia,iL,is,ic,5);

                if isfinite(py_err) && isfinite(oct_err)
                    max_derr = max(max_derr, abs(py_err - oct_err));
                end
                max_dspur = max(max_dspur, abs(py_spur - oct_spur));
                max_dres  = max(max_dres,  abs(py_res  - oct_res));

                if isinf(C.snr_list(is)), lbl = 'inf';
                else, lbl = sprintf('%g', C.snr_list(is)); end
                fprintf('     %-10s  %4d  %5s   %8.2f deg  %8.2f deg  %8.1f%%  %8.1f%%\n', ...
                        C.arms{ia}, Lf, lbl, py_err, oct_err, py_res, oct_res);
            end
        end
    end
end

fprintf('\n');
% The two finders are the same algorithm on the same grid, so the counts must
% match exactly. The angles get a small tolerance only because a symmetric
% crossing's two lobes are a near-tie that the two languages can break
% differently -- see helpers/fODF_peak_score.m.
nfail = nfail + ~check('correct-count fraction reproduces', max_dres < 1e-9, ...
    sprintf('max diff = %.2e %%', max_dres));
nfail = nfail + ~check('spurious-peak count reproduces', max_dspur < 1e-9, ...
    sprintf('max diff = %.2e', max_dspur));
nfail = nfail + ~check('mean angular error reproduces', max_derr < 1e-6, ...
    sprintf('max diff = %.2e deg', max_derr));

fprintf('\n');
if nfail == 0
    fprintf('All checks passed. The .m file can use C.scores directly, and\n');
    fprintf('C.signal to put SMI on the same realisations as the CSD arms.\n');
else
    fprintf('%d check(s) FAILED.\n', nfail);
end
end

% =====================================================================
function ok = check(label, cond, detail)
VERDICT = {'** FAILED **', 'ok'};
ok = logical(cond);
fprintf('   CHECK %-52s %s', label, VERDICT{1+ok});
if nargin > 2 && ~isempty(detail), fprintf('   %s', detail); end
fprintf('\n');
end
