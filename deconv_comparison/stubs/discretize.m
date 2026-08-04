function idx = discretize(x, edges)
% Minimal stand-in for MATLAB's discretize, which Octave does not provide.
% Used by SMI.StandardModel_MLfit_RotInvs. Values outside the edges give NaN,
% the last bin is closed on the right, as MATLAB does.
edges = edges(:)';
idx = nan(size(x));
for k = 1:numel(edges)-1
    if k < numel(edges)-1
        m = x >= edges(k) & x < edges(k+1);
    else
        m = x >= edges(k) & x <= edges(k+1);
    end
    idx(m) = k;
end
end
