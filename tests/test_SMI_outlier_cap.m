% test_SMI_outlier_cap
%
% Integration test: the outlier cap as a flag on SMI.fit. Self contained --
% builds its own protocol and synthetic data, no files needed. Runs in MATLAB,
% and in Octave provided the three usual shims are on the path -- round(x,n),
% discretize(x,edges) and datetime(), which SMI.fit needs and Octave lacks
% (see "README for Claude.md" section 4).
%
% Checks the invariants this codebase holds itself to:
%   1. out.plm / out.pl / out.kernel are bit identical with the flag on or off
%   2. the flag off path is bit identical to a fit that never sets the option
%   3. planted spikes are capped and orientation is preserved exactly
%   4. a contiguous bright region is NOT capped
%   5. the cap and the modulation compose in the documented order
clear; clc
test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
addpath(repo_root);
addpath(fullfile(repo_root, 'helpers'));
if exist('OCTAVE_VERSION','builtin')
    addpath(fullfile(repo_root, 'deconv_comparison', 'stubs'));
end
if exist('OCTAVE_VERSION','builtin'), pkg load statistics; pkg load image; end
TF = {'FAIL','PASS'};

LMAX = 6; CS = 1; D_FW = 3; SNR = 30;
H  = fODF_modulation_helpers();
dq = H.dirs(2000);

% ---- protocol, generated here so the test needs no data
shells = [0 1 2 3];  ndir = [9 60 60 60];
bvals = []; bvecs = [];
for k = 1:numel(shells)
    if shells(k) == 0
        d = zeros(ndir(k),3);
    else
        d = SMI.GetUniformHemisphereDirs(ndir(k));
    end
    bvals = [bvals shells(k)*ones(1,ndir(k))]; %#ok<AGROW>
    bvecs = [bvecs; d];                        %#ok<AGROW>
end
Ndwi = numel(bvals);

% ---- a smoothly varying single-fibre field
sz = [8 8 4]; NV = prod(sz);
kern = [0.60 2.0 2.0 0.50 0.05];
S = zeros(NV, Ndwi);
for iv = 1:NV
    [ix,iy,iz] = ind2sub(sz, iv);
    th = 0.35*ix + 0.2*iy; ph = 0.25*iz;
    n  = [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
    p  = H.watson_plm(dq, n, 16, LMAX, CS);
    S(iv,:) = H.signal(p, [kern 1 1], bvals, ones(1,Ndwi), zeros(1,Ndwi), ...
                       bvecs, LMAX, CS, D_FW)';
end
sigma = 1/SNR;
rand('seed',99); randn('seed',99);
dwi = reshape(sqrt((S + sigma*randn(size(S))).^2 + (sigma*randn(size(S))).^2), [sz Ndwi]);

opts = struct();
opts.b = bvals; opts.dirs = bvecs; opts.sigma = sigma*ones(sz);
opts.mask = true(sz); opts.compartments = {'IAS','EAS','FW'};
opts.NoiseBias = 'Rician'; opts.Lmax = [0 LMAX LMAX LMAX];
opts.CS_phase = CS; opts.D_FW = D_FW; opts.flag_fit_fODF = 1;
opts.fODF_regularization = struct('flag_nonneg',1,'lambda_nonneg',10,'lambda_tikhonov',0.3);

out_off  = SMI.fit(dwi, opts);                        % option never set
o2 = opts; o2.fODF_outlier = struct('flag_cap',0);
out_off2 = SMI.fit(dwi, o2);                          % flag explicitly off

fprintf('\n=== SMI.fit outlier cap integration ===\n\n');
fprintf('[%s] flag off is bit identical to never setting the option (plm)\n', ...
        TF{isequaln(out_off.plm,out_off2.plm)+1});
fprintf('[%s]   ... and pl, and kernel\n', ...
        TF{(isequaln(out_off.pl,out_off2.pl) && isequaln(out_off.kernel,out_off2.kernel))+1});
fprintf('[%s] flag off leaves out.plm_capped absent\n', ...
        TF{(~isfield(out_off2,'plm_capped'))+1});

% ---- plant spikes into a fitted result and cap it post hoc
spikes = [3 3 2; 6 5 3; 4 6 2];
facs   = [40 5e3 1e11];
outS = out_off;
for k = 1:size(spikes,1)
    outS.plm(spikes(k,1),spikes(k,2),spikes(k,3),:) = ...
        outS.plm(spikes(k,1),spikes(k,2),spikes(k,3),:) * facs(k);
end
bx = 6:8; by = 6:8; bz = 2:3;                          % contiguous bright block
outS.plm(bx,by,bz,:) = outS.plm(bx,by,bz,:) * 1.06;
outS.mask = opts.mask; outS.CS_phase = CS;

[plm_cap, ci] = SMI.cap_fODF_outliers(outS, struct('flag_cap',1));
fprintf('\ncapped %d of %d (%.2f%%)  [rel %d, abs %d]  peak max %.4g -> %.4g\n', ...
        ci.Ncap, NV, 100*ci.fraction, ci.Nrel, ci.Nabs, ci.peak_max_before, ci.peak_max_after);

ok = true;
for k = 1:size(spikes,1)
    if ~ci.flagged(spikes(k,1),spikes(k,2),spikes(k,3)), ok = false; end
end
fprintf('[%s] all %d planted spikes flagged\n', TF{ok+1}, size(spikes,1));
blk = ci.flagged(bx,by,bz);
fprintf('[%s] contiguous bright block NOT flagged (%d of %d)\n', ...
        TF{(~any(blk(:)))+1}, nnz(blk), numel(blk));

d = abs(double(plm_cap) - double(outS.plm)); d(~isfinite(d)) = 0;
chg = squeeze(any(d > 0, 4));
fprintf('[%s] only flagged voxels changed (%d changed, %d flagged)\n', ...
        TF{isequal(chg,ci.flagged)+1}, nnz(chg), nnz(ci.flagged));

dirs = SMI.GetUniformHemisphereDirs(500);
Y = SMI.get_even_SH(dirs,LMAX,CS); Y = Y(:,2:end);
L_all = repelem(0:2:LMAX,2*(0:2:LMAX)+1); nrm = sqrt((2*L_all(2:end)+1)/(4*pi));
% Assert the peak lands on the SAME sphere vertex. That is exact; comparing
% angles is not, because acosd cannot resolve below ~1e-5 deg near a dot
% product of 1 and returns floating point noise there.
worst = 0; same_vertex = true;
for k = 1:size(spikes,1)
    p = spikes(k,:);
    a0 = Y*(squeeze(outS.plm(p(1),p(2),p(3),:))'.*nrm)';
    a1 = Y*(squeeze(plm_cap (p(1),p(2),p(3),:))'.*nrm)';
    [~,i0] = max(a0); [~,i1] = max(a1);
    if i0 ~= i1, same_vertex = false; end
    worst = max(worst, acosd(min(abs(dirs(i0,:)*dirs(i1,:)'),1)));
end
fprintf('[%s] peak lands on the identical sphere vertex (angle %.3g deg)\n', ...
        TF{same_vertex+1}, worst);
fprintf('[%s] nothing left above the ceiling 1.0 (max after %.4f)\n', ...
        TF{(ci.peak_max_after <= 1+1e-9)+1}, ci.peak_max_after);

% ---- flag on through SMI.fit, and cap + modulation together. The ceiling is
% deliberately set low here so the cap is forced to act on ordinary voxels.
o4 = opts; o4.fODF_outlier = struct('flag_cap',1,'ceiling',0.30);
out_on = SMI.fit(dwi, o4);
fprintf('\n[%s] flag on leaves out.plm identical to flag off\n', ...
        TF{isequaln(out_on.plm,out_off.plm)+1});
fprintf('[%s] flag on adds out.plm_capped and out.fODF_outlier\n', ...
        TF{(isfield(out_on,'plm_capped') && isfield(out_on,'fODF_outlier'))+1});
fprintf('     ceiling 0.30 forced %d caps of %d voxels\n', out_on.fODF_outlier.Ncap, NV);

o5 = o4; o5.fODF_modulation = struct('flag_modulate',1);
out_both = SMI.fit(dwi, o5);
tmp = out_both; tmp.plm = out_both.plm_capped;    % exactly what SMI.fit does
sh_expect = SMI.modulate_fODF(tmp, o5.fODF_modulation);
fprintf('[%s] modulation was applied to the CAPPED fODF (max|d| = %.3g)\n', ...
        TF{(max(abs(sh_expect(:)-out_both.fODF_modulated(:))) < 1e-12)+1}, ...
        max(abs(sh_expect(:)-out_both.fODF_modulated(:))));
