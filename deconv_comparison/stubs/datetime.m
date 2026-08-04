function s = datetime(varargin)
% Stand-in for MATLAB's datetime, used only by SMI.fit's log writer.
s = datestr(now);
end
