% AER E 351 Homework 02 Problem 2.3 Part b
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Solve for M
M = pi / 2; % [rad]

%% Given
e = 0.5; % []
E_0 = M; % [rad]
tol = 10e-4; % [rad]
n = 5; % []

%% Calculations
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