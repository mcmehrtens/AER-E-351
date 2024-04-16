% Spring 2024 AER E 351 Homework 05 Problem 3.)6.3)
% Matthew Mehrtens
clear; clc; close all;

%% Given
Z = 5; % []
K = 2; % []
I_sp = 300; % [s]
g = 9.81; % [m/s^2]

%% Calculations
t_b = I_sp / K * (1 - 1 / Z); % [s]
v_b = I_sp * g * (log(Z) - 1 / K * (1 - 1 / Z)); % [m/s]

%% Output
fprintf( ...
    "Problem 3.)6.3.)\n" + ...
    "t_b = %g s\n" + ...
    "v_b = %g m/s\n", ...
    t_b, v_b);