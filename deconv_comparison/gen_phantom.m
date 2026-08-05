function gen_phantom(SNR, tag)
% gen_phantom(SNR, tag)
%
% A synthetic brain-like voxel population, whose only purpose is to give
% `dwi2response` something to estimate a response FROM.
%
% READ THIS BEFORE QUOTING ANY NUMBER THAT DEPENDS ON IT. The intended design
% (see Reports/REPORT_SMI_deconvolution_MonteCarlo.md section 2) estimates the
% responses from real data with `dwi2response dhollander`. No real data exists
% in this repository, so the responses are estimated from this phantom instead.
% That makes the response functions REALISTIC -- blunted by fibre dispersion,
% contaminated by whatever the selector lets in, carrying the sampling noise of
% a finite voxel selection -- but not REAL. Nothing here tests how the
% selection behaves on real tissue.
%
% The population deliberately contains what response selection has to survive:
% dispersion spread across the white matter, crossings that a single-fibre
% criterion must reject, partial volume at both the grey matter and the
% ventricle boundary, and grey matter that is not quite isotropic.
%
% Signals come from SMI's own forward model, so the phantom and the Monte
% Carlo simulation cannot end up on different physics. The phantom is also
% written as an MRtrix image, since `dwi2response` reads it directly.

more off
IO = binio();
MR = mrtrix_io();
bvals = IO.load('bvals'); bvals = bvals(:)';
bvecs = IO.load('bvecs');
Ndwi  = numel(bvals);
LMAX_GT = 8; CS = 0; D_FW = 3;

H  = fODF_modulation_helpers();
dq = H.dirs(2000);

GRID = [24 24 16];
N    = prod(GRID);
rand('seed',20260803); randn('seed',20260803);

% class table: {share, kind, [f Da Depar Deperp fw] centre, kappa range}
P = { 0.32, 'single', [0.60 2.00 2.00 0.50 0.02], [10 26]
      0.13, 'cross2', [0.60 2.00 2.00 0.50 0.02], [10 26]
      0.04, 'cross3', [0.60 2.00 2.00 0.50 0.02], [10 26]
      0.08, 'single', [0.42 1.90 1.90 0.65 0.12], [ 6 14]   % WM/GM partial vol
      0.06, 'single', [0.35 2.00 2.20 0.70 0.35], [ 8 18]   % WM/CSF partial vol
      0.27, 'gm',     [0.15 1.20 1.00 0.80 0.05], [0.5 1.5]
      0.10, 'csf',    [0.02 2.00 3.00 3.00 0.95], [0 0] };
KINDS = P(:,2);

counts = round(cell2mat(P(:,1))*N);
counts(end) = N - sum(counts(1:end-1));

L_all = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1);
N_l   = sqrt((2*L_all+1)*(4*pi));
Y_dq  = SMI.get_even_SH(dq, LMAX_GT, CS);
Y_b   = SMI.get_even_SH(bvecs, LMAX_GT, CS);
wq    = 4*pi/size(dq,1);
scale = sqrt((2*L_all'+1)/(4*pi));

% ---- pass 1: microstructure and fODF of every voxel -----------------------
kern  = zeros(N, 5);
A     = zeros(N, numel(L_all));      % normalized plm times N_l
label = zeros(N, 1);
iv = 0;
for ip = 1:size(P,1)
    for k = 1:counts(ip)
        iv = iv + 1;
        kv = P{ip,3} .* (1 + 0.10*(2*rand(1,5)-1));      % 10% jitter per voxel
        kv(1) = min(max(kv(1),0),0.95);
        kv(5) = min(max(kv(5),0),0.98);
        if kv(1)+kv(5) > 0.99, kv(5) = 0.99 - kv(1); end
        kap = P{ip,4}(1) + diff(P{ip,4})*rand();
        switch KINDS{ip}
            case 'csf'
                fod = ones(size(dq,1),1);
            case 'gm'
                fod = H.watson(dq, rand_dir(), kap);
            case 'cross2'
                a = rand_dir(); b2 = rotate_about(a, 40 + 50*rand());
                fod = H.watson(dq,a,kap) + H.watson(dq,b2,kap);
            case 'cross3'
                a = rand_dir(); b2 = rotate_about(a, 50 + 40*rand());
                c = rotate_about(a, 50 + 40*rand());
                fod = H.watson(dq,a,kap) + H.watson(dq,b2,kap) + H.watson(dq,c,kap);
            otherwise
                fod = H.watson(dq, rand_dir(), kap);
        end
        fod = fod(:)/(sum(fod(:))*wq);
        plm_all = ((Y_dq'*fod)*wq)./scale;
        A(iv,:)   = (plm_all/plm_all(1))' .* N_l;
        kern(iv,:) = kv;
        label(iv)  = ip;
    end
end

% ---- pass 2: signals ------------------------------------------------------
% One batched call to the kernel per order, rather than one per voxel: the
% quadrature arrays inside RotInv_Kell_wFW_b_beta_TE_numerical are allocated
% once per call, and per voxel that allocation dominates everything else.
S = zeros(N, Ndwi);
for il = 0:2:LMAX_GT
    Kl  = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(il, bvals, ones(1,Ndwi), ...
              zeros(1,Ndwi), [kern ones(N,2)], D_FW);          % [N x Ndwi]
    idx = find(L_all == il);
    S   = S + Kl.*(A(:,idx)*Y_b(:,idx)');
end

sigma = 1/SNR;
Sn = reshape(sqrt((S + sigma*randn(size(S))).^2 + (sigma*randn(size(S))).^2), ...
             [GRID Ndwi]);
IO.save(['phantom_' tag], Sn);
IO.save(['phantom_label_' tag], label);
IO.save(['phantom_kernel_' tag], kern);

mdir = fullfile(fileparts(mfilename('fullpath')), 'mrtrix');
if ~exist(mdir,'dir'), mkdir(mdir); end
MR.write(fullfile(mdir,['phantom_' tag]), Sn, ...
         struct('grad',[bvecs bvals(:)*1000]));
MR.write(fullfile(mdir,'phantom_mask'), ones(GRID), struct('datatype','UInt8'));

fid = fopen(fullfile(IO.dir(),['phantom_classes_' tag '.txt']),'w');
for ip = 1:size(P,1)
    fprintf(fid,'%s f=%.2f fw=%.2f\n', KINDS{ip}, P{ip,3}(1), P{ip,3}(5));
end
fclose(fid);
fprintf('phantom %s: %d voxels at SNR %g, class counts %s\n', ...
        tag, N, SNR, mat2str(counts'));
end

% =====================================================================
function d = rand_dir()
d = randn(1,3); d = d/norm(d);
end

% =====================================================================
function m = rotate_about(n, deg)
n = n(:)'/norm(n);
t = randn(1,3); e = t - (t*n')*n; e = e/norm(e);
m = cosd(deg)*n + sind(deg)*e; m = m/norm(m);
end
