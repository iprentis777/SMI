function H = binio()
% H = binio()
%
% Flat float64 binary exchange with Python. Column-major, with a sidecar
% '<name>.shape' text file. H.save(name,A) / H.load(name).
H = struct();
H.save = @savebin;
H.load = @loadbin;
H.dir  = @datadir;
end

function d = datadir()
d = fullfile(fileparts(mfilename('fullpath')),'data');
if ~exist(d,'dir'), mkdir(d); end
end

function savebin(name,A)
d = datadir();
A = double(A);
fid = fopen(fullfile(d,[name '.bin']),'wb');
fwrite(fid,A,'float64');            % column-major
fclose(fid);
fid = fopen(fullfile(d,[name '.shape']),'w');
fprintf(fid,'%d ',size(A));
fclose(fid);
end

function A = loadbin(name)
d = datadir();
fid = fopen(fullfile(d,[name '.shape']),'r');
shape = fscanf(fid,'%d');
fclose(fid);
fid = fopen(fullfile(d,[name '.bin']),'rb');
A = fread(fid,prod(shape),'float64');
fclose(fid);
shape = shape(:)';
if numel(shape) < 2, shape = [shape 1]; end   % numpy 1-D -> Octave column
A = reshape(A,shape);
end
