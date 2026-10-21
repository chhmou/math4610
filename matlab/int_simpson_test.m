% int_simpson_test.m
% Test composite Simpson's rule on  int_0^5 x e^{-x} dx
% and estimate the order of accuracy p from
%     p ~ log(err(k)/err(k-1)) / log(h(k)/h(k-1))
% Math 4610, Lecture 20

clear; clc;

f     = @(x) x .* exp(-x);
a     = 0;
b     = 5;
exact = 1 - 6*exp(-5);          % int_0^5 x e^{-x} dx

nlist = 2.^(1:10);              % number of subintervals (must be even)
h     = zeros(size(nlist));
errS  = zeros(size(nlist));
errT  = zeros(size(nlist));

for k = 1:length(nlist)
    n       = nlist(k);
    h(k)    = (b - a)/n;
    errS(k) = abs(exact - int_simpson(f, a, b, n));
    errT(k) = abs(exact - int_trapezoid(f, a, b, n));
end

% observed order of accuracy
pS = [NaN, log(errS(2:end)./errS(1:end-1)) ./ log(h(2:end)./h(1:end-1))];
pT = [NaN, log(errT(2:end)./errT(1:end-1)) ./ log(h(2:end)./h(1:end-1))];

fprintf('%6s %12s %14s %8s %14s %8s\n', 'n', 'h', 'err Simpson', 'p', 'err Trap', 'p');
for k = 1:length(nlist)
    fprintf('%6d %12.4e %14.6e %8.3f %14.6e %8.3f\n', ...
            nlist(k), h(k), errS(k), pS(k), errT(k), pT(k));
end

% log-log plot of error vs h
figure;
loglog(h, errS, 'o-', 'LineWidth', 1.5); hold on;
loglog(h, errT, 's-', 'LineWidth', 1.5);
loglog(h, errS(end)*(h/h(end)).^4, 'k--');
loglog(h, errT(end)*(h/h(end)).^2, 'k:');
grid on;
xlabel('h'); ylabel('|error|');
legend('composite Simpson', 'composite Trapezoid', 'O(h^4)', 'O(h^2)', ...
       'Location', 'southeast');
title('\int_0^5 x e^{-x} dx');

% -------------------------------------------------------------------
function S = int_simpson(f, a, b, n)
% composite Simpson's rule with n (even) uniform subintervals
%   h/3 [ f(x0) + 4 f(x1) + 2 f(x2) + 4 f(x3) + ... + 4 f(x_{n-1}) + f(xn) ]
    if mod(n, 2) ~= 0
        error('n must be even for Simpson''s rule');
    end
    h = (b - a)/n;
    x = a + (0:n)*h;
    y = f(x);
    S = h/3 * ( y(1) + y(end) ...
              + 4*sum(y(2:2:end-1)) ...
              + 2*sum(y(3:2:end-2)) );
end

function T = int_trapezoid(f, a, b, n)
% composite Trapezoid rule with n uniform subintervals
%   h/2 [ f(x0) + 2 f(x1) + ... + 2 f(x_{n-1}) + f(xn) ]
    h = (b - a)/n;
    x = a + (0:n)*h;
    y = f(x);
    T = h/2 * ( y(1) + 2*sum(y(2:end-1)) + y(end) );
end
