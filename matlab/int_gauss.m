function I = int_gauss(f, a, b, n, m)
% INT_GAUSS  Gauss-Legendre quadrature of f over [a,b].
%   I = int_gauss(f, a, b, n)     one panel, n+1 nodes
%   I = int_gauss(f, a, b, n, m)  m equal panels, n+1 nodes per panel
%   f must accept a vector argument.
%
%   Each panel [p,q] is mapped to [-1,1] with
%       x = (q-p)/2 * t + (q+p)/2,   dx = (q-p)/2 dt
    if nargin < 5, m = 1; end
    [t, w] = gauss_nodes(n);          % lookup table, not recomputed
    edges = linspace(a, b, m+1);
    I = 0;
    for k = 1:m
        p = edges(k);  q = edges(k+1);
        x = (q-p)/2 * t + (q+p)/2;    % map nodes to [p,q]
        I = I + (q-p)/2 * sum(w .* f(x));
    end
end
