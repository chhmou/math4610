% INT_GAUSS_TEST  Integrate int_0^5 x e^{-x} dx with Gauss quadrature.
%   (1) 1 subinterval, increasing number of nodes
%   (2) 3 nodes, increasing number of intervals (panels)
clear; clc; close all;

f = @(x) x .* exp(-x);
a = 0;  b = 5;
exact = 1 - 6*exp(-5);

% (1) one subinterval, increasing number of nodes (n+1 nodes)
fprintf('(1) 1 subinterval, increasing nodes\n');
fprintf('%8s %8s %16s %12s\n', 'nodes', 'evals', 'I', 'error');
for n = 0:5
    I = int_gauss(f, a, b, n, 1);
    fprintf('%8d %8d %16.12f %12.3e\n', n+1, n+1, I, abs(I-exact));
end

% (2) 3 nodes (n = 2), increasing number of panels
fprintf('\n(2) 3 nodes, increasing panels\n');
fprintf('%8s %8s %16s %12s\n', 'panels', 'evals', 'I', 'error');
panels = 2.^(0:6);
err2 = zeros(size(panels));
for k = 1:numel(panels)
    I = int_gauss(f, a, b, 2, panels(k));
    err2(k) = abs(I - exact);
    fprintf('%8d %8d %16.12f %12.3e\n', panels(k), 3*panels(k), I, err2(k));
end
