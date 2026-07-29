function s = datetime(varargin)
% Shim: only ever used by SMI's log writer, which prints it.
s = datestr(now);
end
