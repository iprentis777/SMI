function y = round(x,n)
% Shim: MATLAB's two-argument round(x,n). Octave's builtin takes one argument.
if nargin < 2
    y = builtin('round',x);
else
    s = 10^n;
    y = builtin('round',x*s)/s;
end
end
