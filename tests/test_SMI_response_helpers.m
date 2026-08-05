% test_SMI_response_helpers
%
% Self-contained tests for helpers/SMI_response_helpers.m, the zonal harmonic view of
% the SMI kernel used by examples/example_SMI_response_shview.m. Needs no data and no
% fit. Runs in MATLAB and in Octave.
%
% Every test checks the helpers against something computed a different way,
% never against a stored number:
%   1  the zonal response reproduces SMI's own forward model for a delta fODF
%   2  r_0 at b = 0 is S0*sqrt(4*pi), the MRtrix response convention
%   3  an MRtrix response file round trips through write/read
%   4  the glyph radius is the absolute amplitude and its colour the signed one
%   5  an isotropic fODF glyph is a sphere of radius 1/(4*pi)
%   6  a sharp fODF glyph peaks along its own axis
%   7  free water moves r_0 and leaves every l >= 2 coefficient alone
%   8  kernel_from_out pulls the requested voxel out of an SMI.fit output
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

H  = SMI_response_helpers();
HF = fODF_modulation_helpers();
CS = 1; D_FW = 3; LMAX = 8; LMAX_F = 6;
kernel = [0.60 2.0 2.0 0.50 0.05];
bvals  = [0 1 2 3];
ok = true(1,8);

% ---- 1: zonal response vs the forward model, delta fODF ------------------
% p_lm of a delta along z is 1 at m = 0 and 0 elsewhere, in the normalized
% convention. In SMI's ordering the m = 0 entry of block l sits at offset l.
plm_delta = zeros(1, LMAX*(LMAX+3)/2);
ofs = 0;
for il = 2:2:LMAX
    plm_delta(ofs + il + 1) = 1;
    ofs = ofs + 2*il + 1;
end
R = H.zh(kernel, bvals, LMAX, D_FW, 1);
dirs = HF.dirs(300);
e1 = 0;
for i = 1:numel(bvals)
    b_i  = bvals(i)*ones(1,size(dirs,1));
    S_fw = HF.signal(plm_delta, [kernel 1 1], b_i, ones(size(b_i)), ...
                     zeros(size(b_i)), dirs, LMAX, CS, D_FW);
    S_zh = H.profile(R(i,:), acos(min(1,max(-1,dirs(:,3)))));
    e1   = max(e1, max(abs(S_fw(:)-S_zh(:))));
end
ok(1) = e1 < 1e-12;
fprintf('1  zonal response vs SMI forward model   max|err| %.2e   %s\n', e1, TF{ok(1)+1});

% ---- 2: the b = 0 row --------------------------------------------------
S0 = 3.7;
R0 = H.zh(kernel, 0, LMAX, D_FW, S0);
e2 = max(abs(R0(1)-S0*sqrt(4*pi)), max(abs(R0(2:end))));
ok(2) = e2 < 1e-12;
fprintf('2  r_0(b=0) = S0*sqrt(4pi), rest zero    max|err| %.2e   %s\n', e2, TF{ok(2)+1});

% ---- 3: response file round trip ---------------------------------------
fname = [tempname '.txt'];
H.write_response(fname, R);
R_back = H.read_response(fname);
delete(fname);
e3 = max(abs(R(:)-R_back(:)));
ok(3) = isequal(size(R_back), size(R)) && e3 < 1e-7;
fprintf('3  MRtrix response file round trip       max|err| %.2e   %s\n', e3, TF{ok(3)+1});

% ---- 4: glyph geometry --------------------------------------------------
[TH,PH] = H.grid(31,41);
amp = H.profile(R(end,:), TH);
[X,Y,Z,C] = H.glyph(amp, TH, PH);
rad = sqrt(X.^2+Y.^2+Z.^2);
e4 = max(max(abs(rad-abs(amp))));
ok(4) = e4 < 1e-12 && max(max(abs(C-amp))) < 1e-12;
fprintf('4  glyph radius = |amplitude|            max|err| %.2e   %s\n', e4, TF{ok(4)+1});

% ---- 5: isotropic fODF is a sphere of radius 1/(4 pi) -------------------
plm0 = zeros(1, LMAX_F*(LMAX_F+3)/2);
[X,Y,Z,~] = H.sh_glyph(plm0, LMAX_F, CS, 21, 31, 1);
e5 = max(max(abs(sqrt(X.^2+Y.^2+Z.^2) - 1/(4*pi))));
ok(5) = e5 < 1e-12;
fprintf('5  isotropic fODF glyph radius 1/(4pi)   max|err| %.2e   %s\n', e5, TF{ok(5)+1});

% ---- 6: a sharp fODF peaks along its own axis ---------------------------
n  = [0.30 -0.50 0.81]; n = n/norm(n);
dq = HF.dirs(2000);
plm_n = HF.mixture_plm(HF.watson(dq, n, 60), dq, LMAX_F, CS);
[~,~,dirs_g] = H.grid(181,241);
[~,~,~,C] = H.sh_glyph(plm_n, LMAX_F, CS, 181, 241, 1);
[~,imax] = max(C(:));
e6 = acosd(min(1,abs(dirs_g(imax,:)*n')));
ok(6) = e6 < 1.5;
fprintf('6  sharp fODF glyph peaks on its axis    off by  %.2f deg   %s\n', e6, TF{ok(6)+1});

% ---- 7: free water is isotropic ----------------------------------------
% A kernel with fw = 0.30 must have the same l >= 2 coefficients as the sum of
% its intra-axonal and extra-axonal parts alone, since free water is isotropic.
% f_extra = 1 - f - fw = 0.10 here (SMI.m:2447).
% Its whole effect on r_0 must be exactly fw*exp(-b*D_FW)*sqrt(4*pi) -- which
% at b = 3 is 1.3e-4, i.e. free water has essentially left the signal, so the
% check is made across shells rather than at the highest one alone.
b7 = [1 2 3]; fw7 = 0.30;
Rb = H.zh([0.60 2.0 2.0 0.50 fw7], b7, LMAX, D_FW, 1);
Rc = 0.60*H.zh([1 2.0 0 0 0],      b7, LMAX, D_FW, 1) + ...
     0.10*H.zh([0 2.0 2.0 0.50 0], b7, LMAX, D_FW, 1);
e7 = max(max(abs(Rb(:,2:end)-Rc(:,2:end))));
e7b = max(abs((Rb(:,1)-Rc(:,1)) - fw7*exp(-b7(:)*D_FW)*sqrt(4*pi)));
ok(7) = e7 < 1e-12 && e7b < 1e-12;
fprintf('7  free water moves r_0 only             max|err| %.2e / %.2e   %s\n', ...
        e7, e7b, TF{ok(7)+1});

% ---- 8: kernel_from_out -------------------------------------------------
out = struct();
out.kernel = zeros(2,3,4,7);
for ii = 1:2
  for jj = 1:3
    for kk = 1:4
      out.kernel(ii,jj,kk,:) = [ii jj kk 0.5 0.05 80 60];
    end
  end
end
k8 = H.kernel_from_out(out, [2 3 4]);
ok(8) = isequal(k8, [2 3 4 0.5 0.05]);
fprintf('8  kernel_from_out([2 3 4]) = %-22s %s\n', mat2str(k8), TF{ok(8)+1});

fprintf('\n%d passed, %d failed\n', sum(ok), sum(~ok));
if any(~ok), error('test_SMI_response_helpers: %d failure(s)', sum(~ok)); end
