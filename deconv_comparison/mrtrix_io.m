function H = mrtrix_io()
% H = mrtrix_io()
%
% Read and write MRtrix image files from MATLAB/Octave, so that the SMI side of
% the comparison and the MRtrix side operate on the same bytes.
%
%   H.write(name, A, opts)   write an MRtrix image (.mih header + .dat data)
%   H.read(name)             read an MRtrix .mif or .mih image
%   H.write_grad(name, b, g) write an MRtrix gradient table (x y z b per row)
%
% opts fields (all optional):
%   datatype  'Float32LE' (default), 'Float64LE', 'UInt8'
%   grad      [N x 4] gradient table, embedded as dw_scheme
%   vox       voxel sizes, default all 1
%
% FORMAT. The .mih variant is used rather than .mif: the header is a text file
% and the data a separate raw file, which removes the byte offset arithmetic
% that .mif needs. Data is always written with strides 1,2,3,4, i.e. plain
% column-major with axis 0 fastest, which is what fwrite produces from a
% MATLAB array. MRtrix may write its own outputs with different strides, so
% H.read honours whatever `layout` the file declares and returns the array in
% MATLAB order regardless. Run `mrconvert -strides 1,2,3,4` if you want to see
% the raw file in that order too.
H = struct();
H.write      = @write_image;
H.read       = @read_image;
H.write_grad = @write_grad;
end

% =====================================================================
function write_image(name, A, opts)
if nargin < 3, opts = struct(); end
if ~isfield(opts,'datatype'), opts.datatype = 'Float32LE'; end
sz = size(A);
while numel(sz) < 3, sz(end+1) = 1; end %#ok<AGROW>
if ~isfield(opts,'vox') || isempty(opts.vox)
    opts.vox = ones(1, numel(sz));
end

[d, base, ~] = fileparts(name);
hdr = fullfile(d, [base '.mih']);
dat = fullfile(d, [base '.dat']);

switch opts.datatype
    case 'Float32LE', prec = 'float32'; A = single(A);
    case 'Float64LE', prec = 'float64'; A = double(A);
    case 'UInt8',     prec = 'uint8';   A = uint8(A);
    otherwise, error('mrtrix_io: unsupported datatype %s', opts.datatype);
end

fid = fopen(dat, 'wb');
if fid < 0, error('mrtrix_io: cannot write %s', dat); end
fwrite(fid, A, prec);                       % column-major, axis 0 fastest
fclose(fid);

fid = fopen(hdr, 'w');
fprintf(fid, 'mrtrix image\n');
fprintf(fid, 'dim: %s\n', strjoin(arrayfun(@(x) sprintf('%d',x), sz, ...
                                           'UniformOutput', false), ','));
fprintf(fid, 'vox: %s\n', strjoin(arrayfun(@(x) sprintf('%g',x), opts.vox, ...
                                           'UniformOutput', false), ','));
fprintf(fid, 'layout: %s\n', strjoin(arrayfun(@(k) sprintf('+%d',k), ...
                                     0:numel(sz)-1, 'UniformOutput', false), ','));
fprintf(fid, 'datatype: %s\n', opts.datatype);
fprintf(fid, 'transform: 1,0,0,0\n');
fprintf(fid, 'transform: 0,1,0,0\n');
fprintf(fid, 'transform: 0,0,1,0\n');
if isfield(opts,'grad') && ~isempty(opts.grad)
    for i = 1:size(opts.grad,1)
        fprintf(fid, 'dw_scheme: %.10g,%.10g,%.10g,%.10g\n', opts.grad(i,:));
    end
end
fprintf(fid, 'file: %s 0\n', [base '.dat']);
fprintf(fid, 'END\n');
fclose(fid);
end

% =====================================================================
function A = read_image(name)
% Read a .mif or .mih image. Honours `layout`, so an MRtrix output written with
% any strides comes back in MATLAB order.
[d, base, ext] = fileparts(name);
if isempty(ext), ext = '.mih'; name = fullfile(d,[base ext]); end
fid = fopen(name, 'r');
if fid < 0, error('mrtrix_io: cannot read %s', name); end

dim = []; layout = []; datatype = ''; datafile = ''; offset = 0;
while true
    ln = fgetl(fid);
    % A .mif ends its header with END; a .mih written by mrconvert has no END
    % at all, the header simply runs to the end of the file.
    if ~ischar(ln), break, end
    ln = strtrim(ln);
    if strcmp(ln,'END'), break, end
    k = strfind(ln, ':');
    if isempty(k), continue, end
    key = strtrim(ln(1:k(1)-1));
    val = strtrim(ln(k(1)+1:end));
    switch key
        case 'dim',      dim = str2double(strsplit(val, ','));
        case 'layout',   layout = strsplit(val, ',');
        case 'datatype', datatype = val;
        case 'file'
            parts = strsplit(val, ' ');
            datafile = parts{1};
            if numel(parts) > 1, offset = str2double(parts{2}); end
    end
end
hdr_end = ftell(fid);
fclose(fid);

if strcmp(datafile, '.')
    datafile = name;                       % .mif: data follows the header
    if offset == 0, offset = hdr_end; end
else
    datafile = fullfile(d, datafile);
end

switch datatype
    case {'Float32LE','Float32'}, prec = 'float32';
    case {'Float64LE','Float64'}, prec = 'float64';
    case {'UInt8','Bit'},         prec = 'uint8';
    case {'Int16LE','Int16'},     prec = 'int16';
    case {'UInt16LE','UInt16'},   prec = 'uint16';
    otherwise, error('mrtrix_io: unsupported datatype %s', datatype);
end

fid = fopen(datafile, 'rb');
fseek(fid, offset, 'bof');
raw = fread(fid, prod(dim), ['*' prec]);
fclose(fid);
raw = double(raw);

% `layout` gives, per axis, the position of that axis in the file ordering:
% the axis whose entry is +0 varies fastest. Read in file order, then permute.
sgn = ones(1, numel(dim));
pos = zeros(1, numel(dim));
for k = 1:numel(layout)
    s = strtrim(layout{k});
    sgn(k) = 1 - 2*(s(1) == '-');
    pos(k) = str2double(s(2:end));
end
[~, file_order] = sort(pos);                % axes from fastest to slowest
A = reshape(raw, dim(file_order));
A = permute(A, arrayfun(@(k) find(file_order == k), 1:numel(dim)));
for k = 1:numel(dim)
    if sgn(k) < 0, A = flip(A, k); end
end
end

% =====================================================================
function write_grad(name, bvals, dirs)
% MRtrix gradient table: one row per volume, "x y z b" in scanner coordinates.
% The image transform written by write_image is the identity, so scanner
% coordinates and image coordinates are the same thing here.
fid = fopen(name, 'w');
if fid < 0, error('mrtrix_io: cannot write %s', name); end
for i = 1:numel(bvals)
    fprintf(fid, '%.10g %.10g %.10g %.10g\n', dirs(i,1), dirs(i,2), dirs(i,3), bvals(i));
end
fclose(fid);
end
