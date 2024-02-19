% AER E 351 Homework 02 Problem 2.3 Part d
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Solve for M
M = pi / 2; % [rad]

%% Given
e = 0.5; % []
tol = 10e-4; % [rad]
n = 5; % []

%% Calculations
u = M + e; % [rad]
E_0 = (M * (1 - sin(u)) + u * sin(M)) / (1 + sin(M) - sin(u)); % [rad]

F = @(E) E - e * sin(E) - M;
Fprime = @(E) 1 - e * cos(E);
Fprimeprime = @(E) e * sin(E);

E_old = -inf; % [rad]
E_new = E_0; % [rad]
i = 0; % []
fprintf("E_%d = %.4f rad\n", i, E_new)
while abs(E_new - E_old) > tol
    E_old = E_new; % []

    if Fprime(E_new) >= 0
        Fprime_sign = 1; % []
    else
        Fprime_sign = -1; % []
    end

    E_new = E_new ...
        - n * F(E_new) ...
        / (Fprime(E_new) ...
        + Fprime_sign * sqrt((n - 1)^2 * Fprime(E_new)^2 ...
        - n * (n - 1) * F(E_new) * Fprimeprime(E_new))); % [rad]

    i = i + 1; % []
    fprintf("E_%d = %.4f rad\n", i, E_new);
end

%% Calculate f
E = E_new;
f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)); % [rad]
fprintf("f = %.3f rad = %.1f°\n", f, rad2deg(f));