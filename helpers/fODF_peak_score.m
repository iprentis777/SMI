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
function [nfound, aerr] = score_block(SH, Ye, axes_true, ctx)
% [nfound, aerr] = score_block(SH, Ye, axes_true, ctx)
%
% SH is [nrow x ncoef] spherical harmonic coefficients, Ye the matching basis on
% ctx.dirs. axes_true is a cell array of true fibre axes.
nrow   = size(SH,1);
ntrue  = numel(axes_true);
nfound = zeros(nrow,1);
aerr   = nan(nrow,1);
% Each voxel's own isotropic term, NOT the constant. See the header.
A = SH*Ye' - SH(:,1)/sqrt(4*pi);
for r = 1:nrow
    P = row_peaks(A(r,:), ctx);
    nfound(r) = size(P,1);
    if isempty(P), continue, end
    d = zeros(1,ntrue);
    for kk = 1:ntrue
        d(kk) = acosd(min(abs(P(1,:)*axes_true{kk}(:)), 1));
    end
    aerr(r) = min(d);
end
end
