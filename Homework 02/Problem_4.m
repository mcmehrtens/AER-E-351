% AER E 351 Homework 02 Problem 4
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Given
r_0 = [-4622 -5576 0]; % [km]
v_0 = [5.821 -3.772 0]; % [km/s]
mu_Earth = 3.986e5; % [km^3/s^2]

%% Part a
r_0_mag = norm(r_0); % [km]
v_0_mag = norm(v_0); % [km/s]

% v^2 = mu * (2 / r - 1 / a)
a = (2 / r_0_mag - v_0_mag^2 / mu_Earth)^(-1); % [km]

h = cross(r_0, v_0); % [km^2/s]
h_mag = norm(h); % [km^2/s]

% h = sqrt(mu * a * (1 - e^2))
e = sqrt(1 - h_mag^2 / (mu_Earth * a)); % []

p = h_mag^2 / mu_Earth; % [km]

f_0 = 2*pi - acos(1 / e * (p / r_0_mag - 1)); % [rad]
f_0_deg = rad2deg(f_0); % [°]

%% Part b
E_0 = 2 * pi + 2 * atan(sqrt((1 - e) / (1 + e)) * tan(f_0 / 2)); % [rad]
E_0_deg = rad2deg(E_0); % [°]

delta_t = 10 * 60; % [s]

M = sqrt(mu_Earth / a^3) * delta_t; % [rad]

tol = 10e-16; % [rad]
n = 5; % []

%% Calculations
% u = M + e; % [rad]
% E_0 = (M * (1 - sin(u)) + u * sin(M)) / (1 + sin(M) - sin(u)); % [rad]

% F = @(E) E - e * sin(E) - M;
F = @(E) E - E_0 - e * (sin(E) - sin(E_0)) - M;
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
end

E = E_new; % [rad]
E_deg = rad2deg(E); % [°]

f = a / r_0_mag * ((cos(E) - e) * cos(E_0) + sin(E) * sin(E_0)); % []
g = sqrt(a^3 / mu_Earth) * (sin(E - E_0) - e * (sin(E) - sin(E_0))); % []

r = f * r_0 + g * v_0; % [km]
r_mag = norm(r); % [km]

fdot = -sqrt(mu_Earth * a) / (r_mag * r_0_mag) * sin(E - E_0); % []
gdot = a / r_mag * (cos(E_0 - E) - e * cos(E)); % []

v = fdot * r_0 + gdot * v_0; % [km/s]

%% Output
fprintf( ...
    "r_0 = %g km\n" + ...
    "v_0 = %g km/s\n" + ...
    "a = %g km\n" + ...
    "h (km^2/s):\n", ...
    r_0_mag, v_0_mag, a);
disp(h);
fprintf( ...
    "h = %g km^2/s\n" + ...
    "e = %g\n" + ...
    "p = %g km\n" + ...
    "f = %g rad\n" + ...
    "f = %g°\n" + ...
    "E_0 = %g rad\n" + ...
    "E_0 = %g°\n" + ...
    "M = %g rad\n" + ...
    "E = %g rad\n" + ...
    "E = %g°\n" + ...
    "f = %g\n" + ...
    "g = %g\n" + ...
    "r (km):\n", ...
    h_mag, e, p, f_0, f_0_deg, E_0, E_0_deg, M, E, E_deg, f, g);
disp(r);
fprintf( ...
    "r = %g km\n" + ...
    "fdot = %g\n" + ...
    "gdot = %g\n" + ...
    "v (km/s):\n", ...
    r_mag, fdot, gdot);
disp(v);
