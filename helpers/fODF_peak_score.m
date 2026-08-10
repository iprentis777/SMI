function H = fODF_peak_score()
% H = fODF_peak_score()
%
% The peak finder and scorer used to compare deconvolution arms, extracted so
% that the Octave side and the Python side can be shown to agree rather than
% merely to look similar.
%
%   H.setup(dirs_e, nbr_deg)          precompute the neighbour sets, once
%   H.peaks(A_row, ctx)               peak directions of one isotropic-subtracted row
%   H.score(SH, Ye, axes_true, ctx)   [nfound, aerr] for a block of fODFs
%
% This is the algorithm inline in Step 7 of
% deconv_comparison/notebooks/smi_manuscript_60deg.m, unchanged:
%
%   * evaluate the fODF on a dense direction set;
%   * keep every direction not smaller than any neighbour within `nbr_deg`;
%   * of those, keep the ones whose ANISOTROPIC amplitude is at least
%     `rel` of the largest;
%   * sort by amplitude, then remove any peak within `nbr_deg` of a stronger one;
%   * the angular error is the angle from the LARGEST peak to the nearest true
%     axis.
%
% ONE THING THAT MATTERS AND IS EASY TO GET WRONG. The isotropic part
% subtracted here is each voxel's OWN l = 0 term, not the constant 1/(4*pi).
% They are the same number for an fODF in SMI's p_00 = 1 convention and
% different for an MRtrix FOD, which is unnormalised and whose l = 0 coefficient
% varies per voxel and carries apparent fibre density. Subtracting the constant
% would silently mis-threshold every MRtrix arm.
%
% A KNOWN WEAKNESS, recorded because it is visible in the MSMT-CSD results. The
% angular error uses the largest peak only. For a SYMMETRIC crossing the two
% lobes have near-equal amplitude, so which one is "largest" can be decided by
% the last few ulps, and the reported error can jump between the two lobes'
% errors without the reconstruction changing meaningfully. Read it alongside the
% peak COUNT and the spurious count, never on its own.
%
% Lives in its own file rather than as a local function because MATLAB requires
% local functions at the END of a script while Octave cannot call them there at
% all. Same reasoning as helpers/fODF_modulation_helpers.m.
H = struct();
H.setup = @setup_ctx;
H.peaks = @row_peaks;
H.score = @score_block;
end

% =====================================================================
function ctx = setup_ctx(dirs_e, nbr_deg, rel)
% ctx = setup_ctx(dirs_e, nbr_deg, rel)
%
% Precompute the neighbour set of every direction, including itself. This
% replaces the max(a .* nbr, [], 1) form, which built a full ND x ND array for
% every voxel and would dominate the runtime at NREP = 1000.
if nargin < 2 || isempty(nbr_deg), nbr_deg = 12;  end
if nargin < 3 || isempty(rel),     rel     = 0.30; end
ctx = struct();
ctx.dirs = dirs_e;
ctx.ND   = size(dirs_e,1);
ctx.cosN = cosd(nbr_deg);
ctx.rel  = rel;
ctx.nbr  = cell(1, ctx.ND);
for j = 1:ctx.ND
    ctx.nbr{j} = find((dirs_e*dirs_e(j,:)') > ctx.cosN);
end
end

% =====================================================================
function P = row_peaks(a, ctx)
% P = row_peaks(a, ctx)
%
% Peak directions of one row of anisotropic amplitudes, strongest first.
a = a(:)';
Amax = zeros(1, ctx.ND);
for j = 1:ctx.ND
    Amax(j) = max(a(ctx.nbr{j}));
end
lm = find((a > 0) & (a >= Amax));
if isempty(lm), P = zeros(0,3); return, end
lm = lm(a(lm) >= ctx.rel*max(a(lm)));
[~, o] = sort(a(lm), 'descend');
lm = lm(o);
P = ctx.dirs(lm,:);
sel = true(size(P,1),1);
for i = 1:size(P,1)
    if ~sel(i), continue, end
    dup = abs(P*P(i,:)') > ctx.cosN;
    dup(i) = false;
    sel(dup) = false;
end
P = P(sel,:);
end

% =====================================================================
function [nfound, aerr, apk, a0] = score_block(SH, Ye, axes_true, ctx)
% [nfound, aerr, apk, a0] = score_block(SH, Ye, axes_true, ctx)
%
% SH is [nrow x ncoef] spherical harmonic coefficients, Ye the matching basis on
% ctx.dirs. axes_true is a cell array of true fibre axes.
%
% apk is the ANISOTROPIC amplitude at the primary peak -- the fODF amplitude
% there minus that voxel's own isotropic term -- and a0 is the isotropic term
% itself. Together they are what a tractography algorithm actually thresholds
% on, and they are the only outputs here that carry the fODF's SCALE.
%
% *Amplitude is comparable WITHIN an arm, not across arms, and that is not a
% limitation of this code but of the quantity.* An SMI fODF has p_00 = 1 and
% integrates to 1 in every voxel; an MRtrix FOD is unnormalised and its
% amplitude carries apparent fibre density. So `apk` for SMI and `apk` for
% MSMT-CSD are different physical quantities and must never be plotted on one
% axis as though they were. What IS comparable, and is the number worth
% reporting, is the RATIO of an arm's amplitude between conditions -- edema
% against healthy, or one SNR against another -- because the arm's own
% convention cancels.
%
% The neighbourhood maximum is taken for the WHOLE BLOCK at once, one pass per
% direction, rather than per row. At NREP = 1000 the per-row form costs
% nrow*ND small gathers and dominates the runtime of the whole simulation; this
% form costs ND gathers of [nrow x k]. The two select identical peaks.
nrow   = size(SH,1);
ntrue  = numel(axes_true);
nfound = zeros(nrow,1);
aerr   = nan(nrow,1);
apk    = nan(nrow,1);
a0     = SH(:,1)/sqrt(4*pi);        % the isotropic term, per voxel

% Each voxel's own isotropic term, NOT the constant. See the header.
A = SH*Ye' - SH(:,1)/sqrt(4*pi);

Amax = zeros(nrow, ctx.ND);
for j = 1:ctx.ND
    Amax(:,j) = max(A(:, ctx.nbr{j}), [], 2);
end
ismax = (A > 0) & (A >= Amax);

for r = 1:nrow
    lm = find(ismax(r,:));
    if isempty(lm), continue, end
    a  = A(r,:);
    lm = lm(a(lm) >= ctx.rel*max(a(lm)));
    [~, o] = sort(a(lm), 'descend');
    lm = lm(o);
    P  = ctx.dirs(lm,:);
    sel = true(size(P,1),1);
    for i = 1:size(P,1)
        if ~sel(i), continue, end
        dup = abs(P*P(i,:)') > ctx.cosN;
        dup(i) = false;
        sel(dup) = false;
    end
    P = P(sel,:);
    nfound(r) = size(P,1);
    if isempty(P), continue, end
    % The anisotropic amplitude AT the primary peak. lm has already been
    % sorted by amplitude and pruned, so lm(sel) in that order gives the
    % surviving peaks strongest first; a(lm(1)) after the same reordering is
    % the value at P(1,:).
    lmk = lm(sel);
    apk(r) = a(lmk(1));
    d = zeros(1,ntrue);
    for kk = 1:ntrue
        d(kk) = acosd(min(abs(P(1,:)*axes_true{kk}(:)), 1));
    end
    aerr(r) = min(d);
end
end
