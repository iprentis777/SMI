function H = SMI_response_helpers()
% H = SMI_response_helpers()
%
% Helpers for viewing the SMI response kernel as zonal harmonics, the way
% MRtrix3's `shview` displays a response function. Used by
% example_SMI_response_shview.m.
%
%   H.Kell(kernel,b,Lmax,D_FW)        rotational invariants K_l(b) of the kernel
%   H.zh(kernel,b,Lmax,D_FW,S0)       zonal harmonic (m=0) coefficients per shell
%   H.profile(r,theta)                response amplitude R(theta) from zonal r_l
%   H.grid(Nth,Nph)                   (theta,phi) mesh and the direction list
%   H.glyph(amp,TH,PH)                shview-style glyph vertices from amplitudes
%   H.zh_glyph(r,Nth,Nph)             glyph of an axially symmetric response
%   H.sh_glyph(plm,Lmax,CS,Nth,Nph)   glyph of a general (fODF) plm vector
%   H.write_response(file,R)          write an MRtrix-format response .txt
%   H.read_response(file)             read one back
%   H.kernel_from_out(out,idx,opts)   pull a [f Da Depar Deperp fw] kernel out of
%                                     an SMI.fit output struct
%
% CONVENTIONS. SMI's forward model (SMI.m:818-820, and
% fODF_modulation_helpers.m) is
%
%     S(u)/S0 = sum_{lm} K_l(b) p_lm Y_lm(u) sqrt((2l+1)*4*pi)
%
% with p_lm the normalized fODF coefficients, p_00 = 1. A single fibre is a
% delta fODF, whose normalized coefficients are p_l0 = 1 for every even l, so
% the response of one fibre pointing along z is
%
%     R(theta) = sum_l K_l(b) (2l+1) P_l(cos theta)
%              = sum_l r_l Y_l0(theta),    r_l = K_l(b) sqrt((2l+1)*4*pi)
%
% and r_l is exactly what MRtrix stores in a response function text file: one
% row per shell, one column per even l, amplitude reconstructed as
% sum_l r_l Y_l0. That is why the file written by H.write_response can be
% handed straight to `shview`.
%
% These live in their own file rather than as local functions at the end of the
% example because MATLAB requires local functions at the END of a script while
% Octave cannot call them there at all. Same reasoning as
% fODF_modulation_helpers.m.
H = struct();
H.Kell           = @kernel_Kell;
H.zh             = @kernel_zh;
H.profile        = @zh_profile;
H.grid           = @sphere_grid;
H.glyph          = @amp_glyph;
H.zh_glyph       = @zh_glyph;
H.sh_glyph       = @sh_glyph;
H.write_response = @write_response;
H.read_response  = @read_response;
H.kernel_from_out = @kernel_from_out;
end

% =====================================================================
function K = kernel_Kell(kernel,b,Lmax,D_FW)
% K = kernel_Kell(kernel,b,Lmax,D_FW)
%
% [Nb x (Lmax/2+1)] rotational invariants of the Standard Model kernel, one
% column per even l. kernel is [f Da Depar Deperp fw] or [.. T2a T2e].
if nargin < 4 || isempty(D_FW), D_FW = 3; end
b = b(:)';
kernel = kernel(:)';
K = zeros(numel(b),Lmax/2+1);
for il = 0:2:Lmax
    K(:,il/2+1) = SMI.RotInv_Kell_wFW_b_beta_TE_numerical( ...
        il, b, ones(size(b)), zeros(size(b)), kernel, D_FW);
end
end

% =====================================================================
function [R,K] = kernel_zh(kernel,b,Lmax,D_FW,S0)
% [R,K] = kernel_zh(kernel,b,Lmax,D_FW,S0)
%
% Zonal harmonic coefficients of the single fibre response, [Nb x (Lmax/2+1)],
% in the same convention as an MRtrix response function file. S0 scales the
% whole thing (default 1, i.e. the response of the b=0 normalized signal).
if nargin < 4 || isempty(D_FW), D_FW = 3; end
if nargin < 5 || isempty(S0),   S0 = 1;   end
K = kernel_Kell(kernel,b,Lmax,D_FW);
l = 0:2:Lmax;
R = S0 * K .* repmat(sqrt((2*l+1)*4*pi), size(K,1), 1);
end

% =====================================================================
function A = zh_profile(r,theta)
% A = zh_profile(r,theta)
%
% Amplitude of an axially symmetric function with zonal coefficients r
% (one per even l, l = 0,2,...) at polar angles theta (radians, from the
% symmetry axis). A has the shape of theta.
r = r(:)';
Lmax = 2*(numel(r)-1);
x = cos(theta);
A = zeros(size(theta));
for il = 0:2:Lmax
    A = A + r(il/2+1)*sqrt((2*il+1)/(4*pi))*legendre_P(il,x);
end
end

% =====================================================================
function P = legendre_P(l,x)
% Legendre polynomial P_l(x), even orders only, matching SMI.m:2470-2480.
switch l
    case 0, P = ones(size(x));
    case 2, P = 1.5*x.^2 - 0.5;
    case 4, P = (35*x.^4 - 30*x.^2 + 3)/8;
    case 6, P = (231*x.^6 - 315*x.^4 + 105*x.^2 - 5)/16;
    case 8, P = (6435*x.^8 - 12012*x.^6 + 6930*x.^4 - 1260*x.^2 + 35)/128;
    otherwise, error('legendre_P: only even l up to 8 are supported')
end
end

% =====================================================================
function [TH,PH,dirs] = sphere_grid(Nth,Nph)
% [TH,PH,dirs] = sphere_grid(Nth,Nph)
%
% Parametric (theta,phi) mesh covering the whole sphere, plus the same
% directions as an [Nth*Nph x 3] unit vector list. The mesh is what surf()
% needs; the list is what SMI.get_even_SH needs.
th = linspace(0,pi,Nth);
ph = linspace(0,2*pi,Nph);
[PH,TH] = meshgrid(ph,th);
dirs = [sin(TH(:)).*cos(PH(:)), sin(TH(:)).*sin(PH(:)), cos(TH(:))];
end

% =====================================================================
function [X,Y,Z,C] = amp_glyph(amp,TH,PH,scale)
% [X,Y,Z,C] = amp_glyph(amp,TH,PH,scale)
%
% shview-style glyph: the radius is the ABSOLUTE amplitude and the colour is
% the SIGNED amplitude, so negative lobes are visible as a colour change
% rather than being folded silently into the surface. `scale` multiplies the
% radius (default 1). amp must have the size of TH.
if nargin < 4 || isempty(scale), scale = 1; end
amp = reshape(amp,size(TH));
R = scale*abs(amp);
X = R.*sin(TH).*cos(PH);
Y = R.*sin(TH).*sin(PH);
Z = R.*cos(TH);
C = amp;
end

% =====================================================================
function [X,Y,Z,C] = zh_glyph(r,Nth,Nph,scale)
% [X,Y,Z,C] = zh_glyph(r,Nth,Nph,scale)
%
% Glyph of an axially symmetric response with zonal coefficients r, i.e. a
% surface of revolution about z. This is what `shview` draws for one row of a
% response function file.
if nargin < 2 || isempty(Nth), Nth = 121; end
if nargin < 3 || isempty(Nph), Nph = 121; end
if nargin < 4 || isempty(scale), scale = 1; end
[TH,PH] = sphere_grid(Nth,Nph);
[X,Y,Z,C] = amp_glyph(zh_profile(r,TH),TH,PH,scale);
end

% =====================================================================
function [X,Y,Z,C] = sh_glyph(plm,Lmax,CS_phase,Nth,Nph,scale)
% [X,Y,Z,C] = sh_glyph(plm,Lmax,CS_phase,Nth,Nph,scale)
%
% Glyph of a general fODF given in SMI's normalized plm convention (l = 2..Lmax,
% the l=0 term is implicit and contributes the fixed 1/(4*pi) floor). Same
% renderer as the response, so a response and an fODF can be put side by side
% and compared without a convention change in between.
if nargin < 4 || isempty(Nth), Nth = 91;  end
if nargin < 5 || isempty(Nph), Nph = 121; end
if nargin < 6 || isempty(scale), scale = 1; end
[TH,PH,dirs] = sphere_grid(Nth,Nph);
Ylm = SMI.get_even_SH(dirs,Lmax,CS_phase);
L   = repelem(0:2:Lmax,2*(0:2:Lmax)+1)';
amp = Ylm*([1; plm(:)].*sqrt((2*L+1)/(4*pi)));
[X,Y,Z,C] = amp_glyph(amp,TH,PH,scale);
end

% =====================================================================
function write_response(file,R)
% write_response(file,R)
%
% Write [Nshell x Ncoeff] zonal harmonic coefficients as an MRtrix response
% function text file: one row per shell in ascending b order, one column per
% even l starting at 0. `shview <file>` renders exactly these rows.
fid = fopen(file,'w');
if fid < 0, error('write_response: cannot open %s',file); end
for i = 1:size(R,1)
    fprintf(fid,'%.8g',R(i,1));
    fprintf(fid,' %.8g',R(i,2:end));
    fprintf(fid,'\n');
end
fclose(fid);
end

% =====================================================================
function R = read_response(file)
% R = read_response(file)
%
% Read an MRtrix response function text file, skipping '#' comment lines.
fid = fopen(file,'r');
if fid < 0, error('read_response: cannot open %s',file); end
R = [];
while true
    ln = fgetl(fid);
    if ~ischar(ln), break, end
    ln = strtrim(ln);
    if isempty(ln) || ln(1) == '#', continue, end
    R = [R; sscanf(ln,'%f')'];  %#ok<AGROW>
end
fclose(fid);
end

% =====================================================================
function k = kernel_from_out(out,idx,flag_T2)
% k = kernel_from_out(out,idx,flag_T2)
%
% Pull one [f Da Depar Deperp fw] kernel out of an SMI.fit output struct.
% `idx` is either a linear index or a [i j k] subscript into out.kernel.
% out.kernel holds [f Da Depar Deperp fw (T2a T2e) (p2 p4 ..)] along dim 4, so
% only the first five maps are needed here; flag_T2 (default 0) appends the two
% T2 maps, which occupy columns 6 and 7 when the fit used multiple TEs.
if nargin < 3 || isempty(flag_T2), flag_T2 = 0; end
sz = size(out.kernel);
K  = reshape(out.kernel,[],sz(4));
if numel(idx) == 3
    idx = sub2ind(sz(1:3),idx(1),idx(2),idx(3));
end
if flag_T2
    k = K(idx,1:7);
else
    k = K(idx,1:5);
end
end
