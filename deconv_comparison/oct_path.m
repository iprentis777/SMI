% oct_path.m -- put SMI, this directory and the Octave shims on the path.
% MATLAB does not need the shims; they are no-ops there because MATLAB already
% provides round(x,n), discretize and datetime.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
addpath(here);
if exist('OCTAVE_VERSION','builtin')
    addpath(fullfile(here,'stubs'));
    pkg load statistics
    warning('off','all');
    more off
end
