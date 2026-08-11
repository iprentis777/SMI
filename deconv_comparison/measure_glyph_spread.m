function R = measure_glyph_spread(nrep)
% R = measure_glyph_spread([nrep])
%
% How much does a drawn fODF glyph's radius vary BETWEEN noise realisations at
% one SNR, compared with how much it varies ACROSS SNR?
%
%   octave-cli --no-gui -q measure_glyph_spread.m
%   >> measure_glyph_spread          % 32 realisations per SNR, ~4 minutes
%
% This is the instrument behind "README for Claude" section 6.6. It exists
% because Figure 7 of smi_manuscript_60deg.m draws ONE arbitrary realisation per
% SNR -- find(snr_id == is & cond_id == ic, 1) -- and its glyphs come out about
% twice the size at SNR 100 and inf as at SNR 20-50, non-monotonically. The
% question that produced this file was whether that is an SNR effect.
%
% For SSST-CSD it is not: the median radius moves 1.52x across the whole sweep
% while the spread WITHIN a single SNR reaches 2.32x, so which panel looks big
% is mostly which noise draw landed first in the array.
%
% For MSMT-CSD in edema it IS -- and much larger than anyone claimed: the median
% radius moves 7.15x, from 0.286 at SNR 5 to 2.043 noise free. That is a real
% result about MSMT-CSD, not a plotting artefact. See section 6.6.
%
% WHAT IT MEASURES, AND WHY THAT QUANTITY
%
% The maximum clamped glyph radius AFTER the p_00 = 1 normalisation Figure 7
% applies -- that is, literally how big the renderer draws the glyph, not a
% proxy for it. Two earlier answers to this question were wrong (section 2.2 of
% the handoff), and both were wrong by measuring something adjacent: the trend
% over SNR of a single fixed seed, and the behaviour of Octave's axis autoscale
% standing in for MATLAB's. Measure the drawn quantity, over the ensemble the
% figure samples from.
%
% All NREP realisations are packed into one image, so it is one dwi2fod call per
% SNR rather than NREP of them.
%
% MRtrix3 must be on the PATH. Nothing here touches the manuscript file.

if nargin < 1 || isempty(nrep), nrep = 32; end
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'oct_path.m'));

MC = mc_config(); H = fODF_modulation_helpers();
RH = SMI_response_helpers(); MR = mrtrix_io();

% the manuscript settings, edema kernel -- the figure the question was about
K_WM  = [0.60 2.0 2.0 0.50 0.02];        % healthy: the RESPONSE both arms use
K_SIM = [0.10 2.4 2.7 1.15 0.35];        % edema:   what is simulated
KAPPA = 16; D_FW = 3; D_GM = 0.8;
LGT = 8; CS = 0; Lf = 6; B0_SNAP = 0.05;
NEG_LAMBDA = 1; NORM_LAMBDA = 1e-3;
AXIS1 = [0.30 -0.50 0.81]; AXIS1 = AXIS1/norm(AXIS1);
ANG = 60;
SNR_LIST = [5 10 20 30 50 100 Inf];
SNR_LAB  = {'5','10','20','30','50','100','inf'};

D = fullfile(here, 'mrtrix');
if ~exist(D, 'dir'), mkdir(D); end

[bvals, bvecs] = MC.load_protocol_file('hcp_real_3shell.txt');
bvals(bvals < B0_SNAP) = 0;
NB = numel(bvals);
dq = H.dirs(3000);

% Watson zonal coefficients, for the dispersion-matched response
plm_w = H.mixture_plm(H.watson(dq, [0 0 1], KAPPA), dq, LGT, CS);
Mv = []; for il = 0:2:LGT, Mv = [Mv, -il:il]; end %#ok<AGROW>
sh_w = [1; plm_w(:)]; PL = sh_w(Mv(:) == 0);

ax2  = MC.rotate_about(AXIS1, ANG);
fodf = H.watson(dq, AXIS1, KAPPA) + H.watson(dq, ax2, KAPPA);
plm  = H.mixture_plm(fodf, dq, LGT, CS);
S    = H.signal(plm, [K_SIM 1 1], bvals, ones(1,NB), zeros(1,NB), bvecs, ...
                LGT, CS, D_FW);

nc = (Lf/2+1)*(Lf+1);
Lg = repelem(0:2:Lf, 2*(0:2:Lf)+1)';
sc = sqrt((2*Lg(2:end)+1)/(4*pi))';

G = [nrep 1 1];
MR.write(fullfile(D,'gsp_mask.mif'), ones(G), struct('datatype','UInt8'));

randn('seed', 20260811);
fprintf('\n=== glyph radius spread: %d realisations per SNR ===\n', nrep);
fprintf('edema kernel, %d deg crossing, Lmax %d, healthy response\n', ANG, Lf);
fprintf('maximum clamped glyph radius, after the p_00 = 1 normalisation\n\n');
fprintf('%-6s | %-32s | %-32s\n', 'SNR', 'SSST-CSD', 'MSMT-CSD');
fprintf('%-6s | %7s %7s %7s %8s | %7s %7s %7s %8s\n', '', ...
        'min', 'med', 'max', 'max/min', 'min', 'med', 'max', 'max/min');
fprintf('%s\n', repmat('-', 1, 78));

R = struct('snr', SNR_LIST, 'rad', {cell(2, numel(SNR_LIST))}, ...
           'arm', {{'SSST-CSD','MSMT-CSD'}});

for is = 1:numel(SNR_LIST)
    snr = SNR_LIST(is);
    Sn = zeros([nrep 1 1 NB]);
    for r = 1:nrep
        if isfinite(snr)
            sg = 1/snr;                    % Rician, sigma = 1/SNR against S0 = 1
            Sn(r,1,1,:) = sqrt((S(:)' + sg*randn(1,NB)).^2 + (sg*randn(1,NB)).^2);
        else
            Sn(r,1,1,:) = S(:)';
        end
    end
    MR.write(fullfile(D,'gsp_dwi.mif'), Sn, struct('grad', [bvecs bvals(:)*1000]));

    % MRtrix's own per-shell b values: this protocol's b values jitter
    [~, txt] = system(sprintf('mrinfo -shell_bvalues "%s/gsp_dwi.mif" 2>/dev/null', D));
    b_mr = sscanf(txt, '%f')';

    r_wm = RH.zh(K_WM, b_mr/1000, Lf, D_FW) .* repmat(PL(1:Lf/2+1)', numel(b_mr), 1);
    RH.write_response(fullfile(D,'gsp_wm.txt'),  r_wm);
    RH.write_response(fullfile(D,'gsp_b3.txt'),  r_wm(end,:));
    RH.write_response(fullfile(D,'gsp_gm.txt'),  exp(-(b_mr(:)/1000)*D_GM)*sqrt(4*pi));
    RH.write_response(fullfile(D,'gsp_csf.txt'), exp(-(b_mr(:)/1000)*D_FW)*sqrt(4*pi));

    chk(sprintf('dwiextract "%s/gsp_dwi.mif" -shells 0,%g "%s/gsp_b3.mif" -force -quiet 2>&1', ...
                D, b_mr(end), D));
    chk(sprintf(['dwi2fod csd "%s/gsp_b3.mif" "%s/gsp_b3.txt" "%s/gsp_ssst.mif" ' ...
                 '-mask "%s/gsp_mask.mif" -force -quiet 2>&1'], D, D, D, D));
    chk(sprintf(['dwi2fod msmt_csd "%s/gsp_dwi.mif" "%s/gsp_wm.txt" "%s/gsp_mw.mif" ' ...
                 '"%s/gsp_gm.txt" "%s/gsp_mg.mif" "%s/gsp_csf.txt" "%s/gsp_mc.mif" ' ...
                 '-mask "%s/gsp_mask.mif" -neg_lambda %g -norm_lambda %g -force -quiet 2>&1'], ...
                D, D, D, D, D, D, D, D, NEG_LAMBDA, NORM_LAMBDA));

    out = zeros(2,4);
    for ja = 1:2
        if ja == 1, V = MR.read(fullfile(D,'gsp_ssst.mif'));
        else,       V = MR.read(fullfile(D,'gsp_mw.mif')); end
        rad = zeros(1,nrep);
        for r = 1:nrep
            m = squeeze(V(r,1,1,1:nc))';
            % EXACTLY what Figure 7 does: onto p_00 = 1, drop the l = 0 term
            plmf = (m(2:end) * ((1/sqrt(4*pi))/m(1))) ./ sc;
            [Xg,Yg,Zg,~] = RH.sh_glyph(plmf, Lf, CS, 48, 48, 1, 'clamp');
            rad(r) = max(sqrt(Xg(:).^2 + Yg(:).^2 + Zg(:).^2));
        end
        R.rad{ja,is} = rad;
        out(ja,:) = [min(rad) median(rad) max(rad) max(rad)/max(min(rad),eps)];
    end
    fprintf('%-6s | %7.3f %7.3f %7.3f %8.2f | %7.3f %7.3f %7.3f %8.2f\n', ...
            SNR_LAB{is}, out(1,1), out(1,2), out(1,3), out(1,4), ...
            out(2,1), out(2,2), out(2,3), out(2,4));
end

fprintf('\n');
for ja = 1:2
    meds   = cellfun(@median, R.rad(ja,:));
    firsts = cellfun(@(v) v(1), R.rad(ja,:));
    within = max(cellfun(@(v) max(v)/min(v), R.rad(ja,:)));
    fprintf(['%s: the MEDIAN radius moves %.2fx across the whole SNR sweep,\n' ...
             '   but within a single SNR the spread reaches %.2fx. The FIRST\n' ...
             '   realisation of each SNR -- what Figure 7 draws -- moves %.2fx.\n\n'], ...
            R.arm{ja}, max(meds)/min(meds), within, max(firsts)/min(firsts));
end
fprintf(['A figure showing one realisation per SNR therefore cannot show an SNR\n' ...
         'effect on glyph size: the within-SNR spread is larger than the trend.\n' ...
         'The tightening of that spread with SNR, on the other hand, IS a real\n' ...
         'and showable effect -- see "README for Claude" section 6.6.\n\n']);
end

% =====================================================================
function chk(cmd)
[st, txt] = system(cmd);
if st ~= 0
    fprintf(2, 'FAILED: %s\n%s\n', cmd, txt);
    error('measure_glyph_spread: an MRtrix call failed. Is mrtrix3 on the PATH?');
end
end
