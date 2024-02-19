% AER E 351 Homework 02 Problem 2.3 Part g
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

E = fsolve(F, E_0, ...
    optimoptions('fsolve', 'OptimalityTolerance', tol, 'Display', 'none'));

fprintf("E = %.3f rad\n", E);