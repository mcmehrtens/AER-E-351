% AER E 351 Homework 01 Problem 01
% Matthew Mehrtens
clear,clc

% Define given variables
r_c_1 = 9000; % [km]
delta_v_a = 400; % [m/s]
delta_v_b = -400; % [m/s]

% Constants
mu_earth = 3.986e5; % [km^3/s^2]
r_earth = 6.37812e3; % [km]

% Calculate circular orbit before instantaneous change in velocity
% v^2 = mu * (2 / r - 1 / a), vis-viva equation
v_c_1 = sqrt(mu_earth / r_c_1); % [km/s]

% Part 1

% Perogee radius is the same as the circular orbit radius after burn
r_p_2 = r_c_1; % [km]

v_p_2 = v_c_1 + delta_v_a / 1000; % [km/s]

% v^2 = mu * (2 / r - 1 / a), vis-viva equation
a_2 = (2 / r_p_2 - v_p_2^2 / mu_earth)^(-1); % [km]

h_2 = r_p_2 * v_p_2; % [km^2/s]

% h = sqrt(mu * a * (1 - e^2))
e_2 = sqrt(1 - h_2^2 / (mu_earth * a_2));

r_a_2 = a_2 * (1 + e_2); % [km]

alt_p_2 = r_p_2 - r_earth; % [km]
alt_a_2 = r_a_2 - r_earth; % [km]

fprintf( ...
    "Part a.)\n" + ...
    "v_c_1 = %g km/s\n" + ...
    "v_p_2 = %g km/s\n" + ...
    "a_2 = %g km\n" + ...
    "h_2 = %g km^2/s\n" + ...
    "e_2 = %g\n" + ...
    "r_a_2 = %g km\n" + ...
    "alt_p_2 = %g km\n" + ...
    "alt_a_2 = %g km\n\n", ...
    v_c_1, v_p_2, a_2, h_2, e_2, r_a_2, alt_p_2, alt_a_2);

% Part 2

% Apogee radius is the same as the circular orbit radius after burn
r_a_2 = r_c_1; % [km]

v_a_2 = v_c_1 + delta_v_b / 1000; % [km/s]

% v^2 = mu * (2 / r - 1 / a), vis-viva equation
a_2 = (2 / r_a_2 - v_a_2^2 / mu_earth)^(-1); % [km]

h_2 = r_a_2 * v_a_2; % [km^2/s]

% h = sqrt(mu * a * (1 - e^2))
e_2 = sqrt(1 - h_2^2 / (mu_earth * a_2));

r_p_2 = a_2 * (1 - e_2); % [km]

alt_p_2 = r_p_2 - r_earth; % [km]
alt_a_2 = r_a_2 - r_earth; % [km]

fprintf( ...
    "Part b.)\n" + ...
    "v_c_1 = %g km/s\n" + ...
    "v_a_2 = %g km/s\n" + ...
    "a_2 = %g km\n" + ...
    "h_2 = %g km^2/s\n" + ...
    "e_2 = %g\n" + ...
    "r_p_2 = %g km\n" + ...
    "alt_p_2 = %g km\n" + ...
    "alt_a_2 = %g km\n", ...
    v_c_1, v_a_2, a_2, h_2, e_2, r_p_2, alt_p_2, alt_a_2);