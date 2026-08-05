% test_fODF_outlier_cap
%
% Tests SMI.cap_fODF_outliers on a spatially structured phantom. Needs no data
% and no fit: it builds plm directly, so it isolates the cap from everything
% else. Runs in MATLAB and in Octave.
%
% Plants isolated spikes of known size, plants a CONTIGUOUS bright block as the
% edema analogue, and checks that
%   1. every planted spike is flagged and brought to the target
%   2. the contiguous block is NOT flagged by the relative test
%   3. untouched voxels come back bit identical
%   4. peak ORIENTATION of a capped voxel is unchanged
%   5. both rules work in isolation
clear; clc
test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
addpath(repo_root);
addpath(fullfile(repo_root, 'helpers'));
if exist('OCTAVE_VERSION','builtin')
    addpath(fullfile(repo_root, 'deconv_comparison', 'stubs'));
end
if exist('OCTAVE_VERSION','builtin'), pkg load statistics; end
TF = {'FAIL','PASS'};

LMAX = 6; CS = 1; NDIRS = 500;
sz   = [12 12 5];
Nl2  = LMAX*(LMAX+3)/2;                 % plm holds l = 2..Lmax only
iso  = 1/(4*pi);
L_all = repelem(0:2:LMAX, 2*(0:2:LMAX)+1);
nrm   = sqrt((2*L_all(2:end)+1)/(4*pi));
dirs  = SMI.GetUniformHemisphereDirs(NDIRS);
Y     = SMI.get_even_SH(dirs,LMAX,CS); Y = Y(:,2:end);

% quadrature grid for building the Watson fODFs
dq  = SMI.GetUniformHemisphereDirs(2000); dq = [dq; -dq];
Ydq = SMI.get_even_SH(dq,LMAX,CS);
wq  = 4*pi/size(dq,1);

% SMI does not recover the ground truth fODF: measured noise free it returns p2
% at 91%, p4 at 62% and p6 at 28% of truth. Applying that damping makes the
% phantom match what the pipeline actually produces (peak ~0.86 rather than the
% ground truth 1.45), which is what the default ceiling was set from. Without it
% every voxel legitimately exceeds a ceiling of 1 -- which is exactly the
% failure mode the first version of this harness exposed.
damp = [0.907*ones(1,5), 0.616*ones(1,9), 0.277*ones(1,13)];

plm_vol = zeros([sz Nl2]);
for ix = 1:sz(1)
  for iy = 1:sz(2)
    for iz = 1:sz(3)
      th = 0.3*ix + 0.2*iy;  ph = 0.25*iz;
      n  = [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
      kap = 16;
      if ix >= 2 && ix <= 5 && iy >= 2 && iy <= 4 && iz >= 2 && iz <= 4
        kap = 0.3;                      % near-isotropic region (CSF analogue)
      end
      w   = exp(kap*(dq*n(:)).^2); w = w/(sum(w)*wq);
      flm = (Ydq'*w)*wq;
      p   = flm(2:end)./sqrt((2*L_all(2:end)'+1)/(4*pi)) / (flm(1)*sqrt(4*pi));
      plm_vol(ix,iy,iz,:) = p(:).*damp(:);
    end
  end
end
mask = true(sz);

% isolated spikes. Scaling plm scales ONLY the anisotropic part, so the planted
% orientation is known to be unchanged and any rotation the cap introduces would
% be its own doing.
spikes  = [8 8 2; 8 5 3; 5 9 4; 3 3 3];
factors = [50 1e4 1e12 150];
for k = 1:size(spikes,1)
  plm_vol(spikes(k,1),spikes(k,2),spikes(k,3),:) = ...
      plm_vol(spikes(k,1),spikes(k,2),spikes(k,3),:) * factors(k);
end

% contiguous bright block: every voxel raised by the same modest factor, peak
% still under the ceiling. Because the neighbours are raised too, the relative
% test must ignore it. This is the edema analogue.
bx = 9:11; by = 6:8; bz = 2:4;
plm_vol(bx,by,bz,:) = plm_vol(bx,by,bz,:) * 1.08;

in  = struct('plm',plm_vol,'mask',mask,'CS_phase',CS);
opt = struct('orders',1.0,'ceiling',1.0,'min_neighbours',6,'Ndirs',NDIRS);
[plm_out, info] = SMI.cap_fODF_outliers(in, opt);

fprintf('\n=== SMI.cap_fODF_outliers ===\n');
fprintf('capped %d of %d voxels (%.2f%%)  [rel %d, abs %d]\n', ...
        info.Ncap, numel(mask), 100*info.fraction, info.Nrel, info.Nabs);
fprintf('peak before: median %.4f  max %.4g\n', info.peak_median_before, info.peak_max_before);
fprintf('peak after : median %.4f  max %.4g\n\n', info.peak_median_after, info.peak_max_after);

ok = true;
for k = 1:size(spikes,1)
  if ~info.flagged(spikes(k,1),spikes(k,2),spikes(k,3)), ok = false; end
end
fprintf('[%s] all %d planted spikes flagged\n', TF{ok+1}, size(spikes,1));

blk = info.flagged(bx,by,bz);
fprintf('[%s] contiguous bright block NOT flagged (%d of %d)\n', ...
        TF{(~any(blk(:)))+1}, nnz(blk), numel(blk));

d   = abs(plm_out - plm_vol);
chg = squeeze(any(d > 0, 4));
fprintf('[%s] only flagged voxels changed (%d changed, %d flagged)\n', ...
        TF{isequal(chg,info.flagged)+1}, nnz(chg), nnz(info.flagged));

fprintf('\n%-12s %13s %12s %12s %9s\n','spike','peak before','peak after','nbr median','ang err');
worst = 0; ok_t = true; same_vertex = true;
for k = 1:size(spikes,1)
  p  = spikes(k,:);
  a0 = Y*(squeeze(plm_vol(p(1),p(2),p(3),:))'.*nrm)';
  a1 = Y*(squeeze(plm_out(p(1),p(2),p(3),:))'.*nrm)';
  [m0,i0] = max(a0); [m1,i1] = max(a1);
  m0 = m0 + iso;     m1 = m1 + iso;
  ang = acosd(min(abs(dirs(i0,:)*dirs(i1,:)'),1));
  worst = max(worst, ang);
  if i0 ~= i1, same_vertex = false; end
  v = [];
  for dx = -1:1
    for dy = -1:1
      for dz = -1:1
        if dx==0 && dy==0 && dz==0, continue, end
        q = p + [dx dy dz];
        if any(q<1) || any(q>sz), continue, end
        v(end+1) = iso + max(Y*(squeeze(plm_vol(q(1),q(2),q(3),:))'.*nrm)'); %#ok<AGROW>
      end
    end
  end
  if abs(m1 - min(opt.ceiling,median(v))) > 1e-4, ok_t = false; end
  fprintf('%-12s %13.4g %12.4f %12.4f %9.3g\n', mat2str(p), m0, m1, median(v), ang);
end
fprintf('\n[%s] capped peak lands on min(neighbourhood median, ceiling)\n', TF{ok_t+1});
% Vertex identity is the exact test; acosd returns floating point noise below
% ~1e-5 deg near a dot product of 1, so an angle tolerance would be flaky.
fprintf('[%s] peak lands on the identical sphere vertex (angle %.3g deg)\n', ...
        TF{same_vertex+1}, worst);
fprintf('[%s] nothing left above the ceiling %.3g (max after %.4f)\n', ...
        TF{(info.peak_max_after <= opt.ceiling+1e-9)+1}, opt.ceiling, info.peak_max_after);

ii = find(~info.flagged); dmax = 0;
for k = 1:numel(ii)
  [a,b,c] = ind2sub(sz, ii(k));
  p0 = iso + max(Y*(squeeze(plm_vol(a,b,c,:))'.*nrm)');
  p1 = iso + max(Y*(squeeze(plm_out(a,b,c,:))'.*nrm)');
  dmax = max(dmax, abs(p0-p1));
end
fprintf('[%s] every unflagged voxel bit identical (max|d| = %.3g over %d)\n', ...
        TF{(dmax==0)+1}, dmax, numel(ii));

% ---- ceiling effectively disabled, to isolate the RELATIVE rule
opt2 = opt; opt2.ceiling = 1e9;
[~, info_r] = SMI.cap_fODF_outliers(in, opt2);
okr = true;
for k = 1:size(spikes,1)
  if ~info_r.flagged(spikes(k,1),spikes(k,2),spikes(k,3)), okr = false; end
end
blkr = info_r.flagged(bx,by,bz);
fprintf('\n--- ceiling disabled (1e9), relative rule alone ---\n');
fprintf('[%s] all %d spikes still caught by the neighbourhood rule alone\n', ...
        TF{okr+1}, size(spikes,1));
fprintf('[%s] contiguous block still not flagged (%d of %d)\n', ...
        TF{(~any(blkr(:)))+1}, nnz(blkr), numel(blkr));
fprintf('[%s] nothing else caught (%d flagged, %d planted)\n', ...
        TF{(info_r.Ncap == size(spikes,1))+1}, info_r.Ncap, size(spikes,1));
