% AER E 351 Homework 03 Problem 3.2a
% Matthew Mehrtens
clear, clc, close all;

%% Given
r = [sqrt(2) / 2
    sqrt(2) / 2
    0];
v = [-sqrt(2) / 2
    sqrt(2) / 2
    0];

%% Calculations
r_mag = norm(r);
v_mag = norm(v);

% a = r / (2 - rv^2 / mu)
a = r_mag / (2 - r_mag * v_mag^2);

% e = [(|v|^2 - mu / |r|) * r - (r . v) * v]
% dot(r, v) % <-- r.v = 0
e = (v_mag^2 - 1 / r_mag) * r;
e_mag = norm(e);

h = cross(r, v);
h_mag = norm(h);

i = acosd(h(3) / h_mag); % [°]

n = cross([0 0 1], h / h_mag);
n_mag = norm(n);

% Since n_y is < 0, we need a correction on omega
Omega = 180 + acosd(n(1) / n_mag);

omega = acosd(dot(n, e) / (n_mag * e_mag)); % [°]

f = acosd(0 / 0); % [°]

%% Output
fprintf( ...
    "Problem 3.2b Solutions:\n" + ...
    "a = %g\n" + ...
    "|e| = %g\n" + ...
    "i = %g°\n" + ...
    "Omega = %g°\n" + ...
    "omega = %g°\n" + ...
    "f = %g°\n", ...
    a, e_mag, i, Omega, omega, f);