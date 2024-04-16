% AER E 351 Homework 03 Problem 1a
% Matthew Mehrtens
clear, clc, close all;

%% Given
R_Earth = 6.37812e3; % [km]
mu_Earth = 3.986e5; % [km^3/s^2]
a = 1.7 * R_Earth; % [km]
e = 0.4; % []
i = deg2rad(20); % [rad]
Omega = deg2rad(30); % [rad]
omega = deg2rad(10); % [rad]
f_0 = deg2rad(60); % [rad]
theta = omega + f_0; % [rad]

%% Calculations
r_mag = a * (1 - e^2) / (1 + e * cos(f_0)); % [km]
r = r_mag * [...
    cos(Omega) * cos(theta) - sin(Omega) * sin(theta) * cos(i)
    sin(Omega) * cos(theta) + cos(Omega) * sin(theta) * cos(i)
    sin(theta) * sin(i)]; % [km]

h = sqrt(mu_Earth * a * (1 - e^2)); % [km^2/s]

v = mu_Earth / h * [...
    -(cos(Omega) * (sin(theta) + e * sin(omega)) + sin(Omega) * (cos(theta) + e * cos(omega)) * cos(i))
    -(sin(Omega) * (sin(theta) + e * sin(omega)) - cos(Omega) * (cos(theta) + e * cos(omega)) * cos(i))
    (cos(theta) + e * cos(omega)) * sin(i)]; % [km/s]

%% Display
fprintf(...
    "|r| = %g km\n" + ...
    "r [km]:\n", ...
    r_mag);
disp(r);
fprintf(...
    "h = %g km^2/s\n" + ...
    "v [km/s]:\n", ...
    h);
disp(v);