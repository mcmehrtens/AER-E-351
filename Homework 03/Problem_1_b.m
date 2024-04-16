% AER E 351 Homework 03 Problem 1b
% Matthew Mehrtens
clear, clc, close all;

%% Given
R_Earth = 6.37812e3; % [km]
mu_Earth = 3.986e5; % [km^3/s^2]
a = 2 * R_Earth; % [km]
e = 0; % []
i = deg2rad(10); % [rad]
Omega = deg2rad(30); % [rad]
theta_0 = deg2rad(20); % [rad]

%% Calculations
r_mag = a; % [km]
r = r_mag * [...
    cos(Omega) * cos(theta_0) - sin(Omega) * sin(theta_0) * cos(i)
    sin(Omega) * cos(theta_0) + cos(Omega) * sin(theta_0) * cos(i)
    sin(theta_0) * sin(i)]; % [km]

v_mag = sqrt(mu_Earth * (2 / r_mag - 1 / a)); % [km/s]
h = sqrt(mu_Earth * a * (1 - e^2)); % [km^2/s]
gamma = acos(h / (r_mag * v_mag)); % [rad]
c = cross(r, [cos(Omega) sin(Omega) 0]); % [km]

syms v_x v_y v_z;
S = solve(...
    v_mag^2 == v_x^2 + v_y^2 + v_z^2, ...
    r_mag * v_mag * sin(gamma) == r(1) * v_x + r(2) * v_y + r(3) * v_z, ...
    0 == c(1) * v_x + c(2) * v_y + c(3) * v_z, ...
    v_x, v_y, v_z);
v = [double(S.v_x(2)) double(S.v_y(2)) double(S.v_z(2))]; % [km/s]

%% Display
fprintf(...
    "|r| = %g km\n" + ...
    "r [km]:\n", ...
    r_mag);
disp(r);
fprintf(...
    "|v| = %g km/s\n" + ...
    "h = %g km^2/s\n" + ...
    "gamma = %g rad\n" + ...
    "c [km]:\n", ...
    v_mag, h, gamma);
disp(c);
fprintf("v [km/s]:\n");
disp(v);