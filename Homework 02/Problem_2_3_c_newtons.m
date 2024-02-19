% AER E 351 Homework 02 Problem 2.3 Part c
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Solve for M
M = pi / 2; % [rad]

%% Given
e = 0.5; % []
tol = 10e-4; % [rad]

%% Calculations
u = M + e; % [rad]
E_0 = (M * (1 - sin(u)) + u * sin(M)) / (1 + sin(M) - sin(u)); % [rad]

F = @(E) E - e * sin(E) - M;
Fprime = @(E) 1 - e * cos(E);

E_old = -inf; % [rad]
E_new = E_0; % [rad]
i = 0; % []
fprintf("E_%d = %.4f rad\n", i, E_new)
while abs(E_new - E_old) > tol
    E_old = E_new; % []
    E_new = E_new - F(E_new) / Fprime(E_new); % [rad]
    i = i + 1; % []
    fprintf("E_%d = %.4f rad\n", i, E_new);
end