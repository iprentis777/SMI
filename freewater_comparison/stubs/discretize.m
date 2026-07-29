function idx = discretize(x,edges)
% Shim: MATLAB's discretize. Returns the index of the bin of edges containing
% each element of x, NaN outside. Last bin is closed on the right.
edges = edges(:)';
idx = nan(size(x));
for k = 1:numel(edges)-1
    if k < numel(edges)-1
        in = x >= edges(k) & x < edges(k+1);
    else
        in = x >= edges(k) & x <= edges(k+1);
    end
    idx(in) = k;
end
end
