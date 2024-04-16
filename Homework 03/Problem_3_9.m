% AER E 351 Homework 03 Problem 3.9
% Matthew Mehrtens
clear, clc, close all;

%% Given
e = 0.5; % []
f = 45; % [°]

%% Problem 3.9a
gamma_a = atand(0.707 / 0.707) - atand((0.707 - 0.293) / 0.707); % [°]
verify_gamma_a = atand(e * sind(f) / (1 + e * cosd(f))); % [°]

%% Problem 3.9b
gamma_min = -asind(e); % [°]

%% Problem 3.9c
p_coeff_c = 1 / (1 + 0.354); % []
verify_p_coeff_c = 1 / (1 + e * cosd(f)); % []

%% Problem 3.9d
p_coeff_d = 1 / (1 - 0.25); % []
verify_p_coeff_d = 1 / (1 + e * cosd(240)); % []

%% Output
fprintf( ...
    "Problem 3.9 Solutions:\n" + ...
    "a:\n" + ...
    "gamma = %g°\n" + ...
    "verification: gamma = %g°\n" + ...
    "b:\n" + ...
    "gamma_min = %g°\n" + ...
    "c:\n" + ...
    "r = %gp\n" + ...
    "verification: r = %gp\n" + ...
    "d:\n" + ...
    "r = %gp\n" + ...
    "verification: r = %gp\n" + ...
    "e:\n" + ...
    "Yes, by using Equation 1.33:\n" + ...
    "\tr = a * (1 - e^2) / (1 + e * cos(f))\t(Eq. 1.33)\n", ...
    gamma_a, verify_gamma_a, gamma_min, p_coeff_c, verify_p_coeff_c, ...
    p_coeff_d, verify_p_coeff_d);