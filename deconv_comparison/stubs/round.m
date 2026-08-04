function y = round(x,n)
% MATLAB's two-argument round(x,n), which Octave does not provide.
% Used by SMI.Group_dwi_in_shells_b_beta_TE.
if nargin < 2, n = 0; end
y = builtin('round', x .* 10.^n) ./ 10.^n;
end
