function C = read_csd_export(data_dir)
% C = read_csd_export()
%
% Read everything notebooks/csd_manuscript_60deg.ipynb exported, so the .m side
% can draw figures over the CSD arms without re-running any of them.
%
% Returns a struct:
%
%   C.preset        'healthy'
%   C.kernel        [f Da Depar Deperp fw] the CSD arms were simulated from
%   C.kappa         Watson concentration of each fibre population
%   C.response_mode 'dispersed' or 'delta' -- which response the binaries got
%   C.protocol      the acquisition file name
%   C.cs_phase      0, i.e. SMI's basis IS MRtrix's
%   C.lmax_gt       ground truth angular order
%   C.angles        crossing angles simulated
%   C.snr_list      the SNR sweep, Inf for the noise-free block
%   C.lmax_list     the angular orders fitted at
%   C.arms          {1 x NARM} cell of arm names, e.g. {'SSST-CSD','MSMT-CSD'}
%   C.sh{ia}{iL}    [NVOX x ncoef] SH coefficients, arm ia at Lmax C.lmax_list(iL)
%   C.scores        [NARM x NLMAX x NSNR x NCOND x 5], last index in the order
%                   C.score_order: pct_correct mean_err std_err median_err spurious
%   C.ceiling       [NARM x NLMAX x NCOND x 2]: peak count and error of the
%                   band-limited truth, i.e. what the angular order allows
%   C.signal        [NVOX x Ndwi] THE NOISY SIGNAL the CSD arms were given
%   C.truth_sh      [NCOND x ncoef] the ground truth at C.lmax_gt
%   C.bvals,C.bvecs the acquisition, as the notebook read it
%   C.axes{ic}      {1 x ntrue} true fibre axes of condition ic
%   C.cond_id       [NVOX x 1] which condition each voxel is, 1-based
%   C.snr_id        [NVOX x 1] which SNR block each voxel is, 1-based
%   C.grid_all      the 3D grid the voxels are laid out on
%
% C.signal is the one that matters most. Fitting SMI to it, rather than to a
% fresh noise draw, is what puts all three arms on THE SAME realisations -- so a
% difference between arms is the method and not the Monte Carlo error of two
% independent noise draws.
%
% Requires the notebook to have been run. Everything lives in data/, which is
% gitignored, so this errors clearly rather than returning stale numbers.

if nargin < 1 || isempty(data_dir)
    data_dir = fullfile(fileparts(mfilename('fullpath')), 'data');
end
man = fullfile(data_dir, 'csd_manifest.txt');
if ~exist(man, 'file')
    error(['read_csd_export: %s not found.\n' ...
           'Run notebooks/csd_manuscript_60deg.ipynb first -- its Step 8 writes it.'], man);
end

IO = binio();
M  = read_manifest(man);

C = struct();
C.preset        = M.preset;
C.kernel        = M.kernel;
C.kappa         = M.kappa;
C.response_mode = M.response_mode;
C.protocol      = M.protocol;
C.b0_snap       = M.b0_snap;
C.cs_phase      = M.cs_phase;
C.lmax_gt       = M.lmax_gt;
C.smoke_test    = M.smoke_test;
C.nrep          = M.nrep;
C.seed          = M.seed;
C.angles        = M.angles;
C.lmax_list     = M.lmax_list;
C.grid_all      = M.grid_all;
C.nvox          = M.nvox;
C.nvox_per_snr  = M.nvox_per_snr;
C.ndwi          = M.ndwi;
C.arms          = M.arms;
C.score_order   = M.score_order;

% snr_list arrives as text so that 'inf' survives the round trip intact.
C.snr_list = zeros(1, numel(M.snr_list_txt));
for i = 1:numel(M.snr_list_txt)
    if strcmpi(M.snr_list_txt{i}, 'inf')
        C.snr_list(i) = Inf;
    else
        C.snr_list(i) = str2double(M.snr_list_txt{i});
    end
end

NCOND = numel(C.angles);
NSNR  = numel(C.snr_list);
NARM  = numel(C.arms);
NL    = numel(C.lmax_list);

% The true fibre axes, per condition.
C.axes = cell(1, NCOND);
for ic = 1:NCOND
    ax = {};
    for j = 0:1
        key = sprintf('axis_cond%d_%d', ic-1, j);
        if isfield(M.raw, key), ax{end+1} = M.raw.(key); end %#ok<AGROW>
    end
    C.axes{ic} = ax;
end

% Voxel layout, rebuilt from the manifest rather than stored, so it cannot
% disagree with the counts.
C.cond_id = repmat(repelem((1:NCOND)', C.nrep, 1), NSNR, 1);
C.snr_id  = repelem((1:NSNR)', C.nvox_per_snr, 1);

C.signal   = IO.load('csd_signal');
C.truth_sh = IO.load('csd_truth_sh');
C.bvals    = IO.load('csd_bvals');
C.bvals    = C.bvals(:)';
C.bvecs    = IO.load('csd_bvecs');
C.scores   = reshape(IO.load('csd_scores'),  [NARM NL NSNR NCOND 5]);
C.ceiling  = reshape(IO.load('csd_ceiling'), [NARM NL NCOND 2]);

C.sh = cell(1, NARM);
for ia = 1:NARM
    C.sh{ia} = cell(1, NL);
    nm = strrep(C.arms{ia}, '-', '');
    for iL = 1:NL
        C.sh{ia}{iL} = IO.load(sprintf('csd_sh_%s_lmax%d', nm, C.lmax_list(iL)));
    end
end

% Consistency, checked here rather than trusted, because a mismatch means the
% manifest and the arrays came from different runs of the notebook.
assert(size(C.signal,1) == C.nvox, ...
       'signal has %d voxels, manifest says %d', size(C.signal,1), C.nvox);
assert(size(C.signal,2) == C.ndwi, ...
       'signal has %d volumes, manifest says %d', size(C.signal,2), C.ndwi);
assert(numel(C.cond_id) == C.nvox, 'voxel layout does not match nvox');
for ia = 1:NARM
    for iL = 1:NL
        assert(size(C.sh{ia}{iL},1) == C.nvox, ...
               '%s Lmax %d has %d voxels, expected %d', C.arms{ia}, ...
               C.lmax_list(iL), size(C.sh{ia}{iL},1), C.nvox);
    end
end

fprintf('read_csd_export: %s kernel, %s response, %d arms (%s)\n', ...
        C.preset, C.response_mode, NARM, strjoin(C.arms, ', '));
fprintf('   %d voxels = %d conditions x %d reps x %d SNR, Lmax %s\n', ...
        C.nvox, NCOND, C.nrep, NSNR, mat2str(C.lmax_list));
if C.smoke_test
    fprintf('   *** SMOKE_TEST run: reduced sweep, indicative numbers only ***\n');
end
end

% =====================================================================
function M = read_manifest(fname)
% Parse the key/value manifest. Everything numeric becomes a row vector;
% snr_list and arms stay as text because 'inf' and 'SSST-CSD' are not numbers.
M = struct(); M.raw = struct();
fid = fopen(fname, 'r');
if fid < 0, error('cannot open %s', fname); end
txt_keys = {'preset','response_mode','protocol'};
while true
    ln = fgetl(fid);
    if ~ischar(ln), break, end
    ln = strtrim(ln);
    if isempty(ln) || ln(1) == '%', continue, end
    sp  = find(ln == ' ', 1);
    if isempty(sp), continue, end
    key = ln(1:sp-1);
    val = strtrim(ln(sp+1:end));
    if any(strcmp(key, txt_keys))
        M.(key) = val;
    elseif strcmp(key, 'snr_list')
        M.snr_list_txt = strsplit(val, ' ');
    elseif strcmp(key, 'arms')
        M.arms = strsplit(val, ' ');
    elseif strcmp(key, 'score_order')
        M.score_order = strsplit(val, ' ');
    else
        M.(key) = sscanf(val, '%f')';
        M.raw.(key) = M.(key);
    end
end
fclose(fid);
end
