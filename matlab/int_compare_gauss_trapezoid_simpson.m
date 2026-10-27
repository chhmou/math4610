% INT_COMPARE_GAUSS_TRAPEZOID_SIMPSON
%   Integrate int_0^5 x e^{-x} dx with composite
%     1. trapezoid   (2 nodes per panel)
%     2. Simpson     (3 nodes per panel)
%     3. Gauss       (3 nodes per panel, n = 2)
%   using the same panels, and compare error vs. function evaluations.
%   (Local functions at the end of a script need MATLAB R2016b or later.)
clear; clc; close all;

f = @(x) x .* exp(-x);
a = 0;  b = 5;
exact = 1 - 6*exp(-5);

panels = 2.^(0:7);
N = numel(panels);
errT = zeros(1,N);  errS = zeros(1,N);  errG = zeros(1,N);
evT  = zeros(1,N);  evS  = zeros(1,N);  evG  = zeros(1,N);

for k = 1:N
    m = panels(k);
    errT(k) = abs(trapezoid_panels(f, a, b, m) - exact);   evT(k) = m + 1;
    errS(k) = abs(simpson_panels(f, a, b, m)   - exact);   evS(k) = 2*m + 1;
    errG(k) = abs(int_gauss(f, a, b, 2, m)     - exact);   evG(k) = 3*m;
end

fprintf('%7s | %6s %11s | %6s %11s | %6s %11s\n', 'panels', ...
        'evalT', 'errT', 'evalS', 'errS', 'evalG', 'errG');
for k = 1:N
    fprintf('%7d | %6d %11.3e | %6d %11.3e | %6d %11.3e\n', panels(k), ...
            evT(k), errT(k), evS(k), errS(k), evG(k), errG(k));
end

figure;
loglog(evT, errT, 'o-', evS, errS, 's-', evG, errG, 'd-', 'LineWidth', 1.5);
grid on;
xlabel('number of function evaluations'); ylabel('error');
legend('Trapezoid', 'Simpson', 'Gauss (3 nodes)', 'Location', 'southwest');
title('\int_0^5 x e^{-x} dx');

% ---------------- local functions ----------------
function I = trapezoid_panels(f, a, b, m)
    x = linspace(a, b, m+1);
    h = (b - a)/m;
    y = f(x);
    I = h * (sum(y) - (y(1) + y(end))/2);
end

function I = simpson_panels(f, a, b, m)
    x = linspace(a, b, 2*m+1);         % panel endpoints and midpoints
    h = (b - a)/m;
    y = f(x);
    I = h/6 * (y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)));
end
