function gen_brain(SNR, variant, tag)
% gen_brain(SNR, variant, tag)
%
% A synthetic "brain" population, purely so that CSD and MSMT-CSD can estimate
% their response functions the way the real tools do -- from the most
% anisotropic voxels in the volume -- instead of being handed the ground truth.
%
% variant 'healthy'  no edema anywhere
%         'edema'    15% of the volume is edematous with the intra-axonal
%                    fraction INTACT (f = 0.60) and half the extra-axonal
%                    water converted to free water (fw = 0.20)
%
% The point of the two variants: free water is essentially fully attenuated by
% b = 3 ms/um^2, so converting hindered extra-axonal water into free water
% RAISES the anisotropy of what survives at high b. An edematous voxel can
% therefore look MORE single-fibre-like than healthy white matter to a
% response estimator that selects on anisotropy.

more off
IO = binio();
bvals = IO.load('bvals')(:)';
bvecs = IO.load('bvecs');
Ndwi  = numel(bvals);
LMAX_GT = 8; CS = 1; D_FW = 3;
H = fODF_modulation_helpers();
dq = H.dirs(2000);

GRID = [20 20 10];
N = prod(GRID);
rand('seed',4242); randn('seed',4242);

% tissue classes: {fraction, kernel, kind}
if strcmp(variant,'edema')
    P = { 0.35, [0.60 2.0 2.0 0.50 0.02], 'single'
          0.20, [0.60 2.0 2.0 0.50 0.02], 'cross'
          0.15, [0.60 2.0 2.0 0.50 0.20], 'single'      % edema, f intact
          0.10, [0.45 2.0 2.0 0.60 0.15], 'single'      % mild partial volume
          0.10, [0.15 1.2 1.0 0.80 0.10], 'gm'
          0.10, [0.02 2.0 3.0 3.00 0.95], 'csf' };
else
    P = { 0.45, [0.60 2.0 2.0 0.50 0.02], 'single'
          0.25, [0.60 2.0 2.0 0.50 0.02], 'cross'
          0.10, [0.45 2.0 2.0 0.60 0.15], 'single'
          0.10, [0.15 1.2 1.0 0.80 0.10], 'gm'
          0.10, [0.02 2.0 3.0 3.00 0.95], 'csf' };
end

counts = round(cell2mat(P(:,1))*N);
counts(end) = N - sum(counts(1:end-1));

% Precompute everything that does not depend on the individual voxel: the SH
% bases, and one kernel per tissue class. Done naively this loop calls
% get_even_SH 8000 times.
L_all = repelem(0:2:LMAX_GT, 2*(0:2:LMAX_GT)+1);
N_l   = sqrt((2*L_all+1)*(4*pi));
Y_dq  = SMI.get_even_SH(dq, LMAX_GT, CS);            % [Ndq x Nlm]
Y_b   = SMI.get_even_SH(bvecs, LMAX_GT, CS);         % [Ndwi x Nlm]
wq    = 4*pi/size(dq,1);
scale = sqrt((2*L_all'+1)/(4*pi));

DES = cell(size(P,1),1);
for ip = 1:size(P,1)
    Kell = zeros(Ndwi, LMAX_GT/2+1);
    for il = 0:2:LMAX_GT
        Kell(:,il/2+1) = SMI.RotInv_Kell_wFW_b_beta_TE_numerical( ...
            il, bvals, ones(1,Ndwi), zeros(1,Ndwi), [P{ip,2} 1 1], D_FW);
    end
    Kmat = Kell(:, repelem(1:(LMAX_GT/2+1), 2*(0:2:LMAX_GT)+1));
    DES{ip} = Kmat.*(Y_b.*N_l);                      % [Ndwi x Nlm]
end

S = zeros(N, Ndwi);
label = zeros(N,1);
iv = 0;
for ip = 1:size(P,1)
    for k = 1:counts(ip)
        iv = iv + 1;
        switch P{ip,3}
            case 'csf'
                fod = ones(size(dq,1),1);
            case 'gm'
                fod = H.watson(dq, rand_dir(), 0.8);
            case 'cross'
                a = rand_dir(); b = rotate_about(a, 50 + 40*rand());
                fod = H.watson(dq,a,12+8*rand()) + H.watson(dq,b,12+8*rand());
            otherwise
                fod = H.watson(dq, rand_dir(), 12 + 8*rand());
        end
        fod = fod(:)/(sum(fod(:))*wq);               % integrates to 1
        flm = (Y_dq'*fod)*wq;
        plm_all = flm./scale;
        S(iv,:) = (DES{ip}*(plm_all/plm_all(1)))';   % p_00 = 1 convention
        label(iv) = ip;
    end
end

sigma = 1/SNR;
Sn = sqrt((S + sigma*randn(size(S))).^2 + (sigma*randn(size(S))).^2);
IO.save(['brain_' tag], reshape(Sn, [GRID Ndwi]));
IO.save(['brain_label_' tag], label);
fprintf('brain %s (%s): %d voxels, classes %s\n', tag, variant, N, mat2str(counts'));
end

function d = rand_dir()
d = randn(1,3); d = d/norm(d);
end

function m = rotate_about(n, deg)
n = n(:)'/norm(n);
t = randn(1,3); e = t - (t*n')*n; e = e/norm(e);
m = cosd(deg)*n + sind(deg)*e; m = m/norm(m);
end
