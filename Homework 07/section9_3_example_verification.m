% Spring 2024 AER E 351 Section 9.3 Example Verification
% Matthew Mehrtens
clear; close all; clc;

%% Equations
emos2kmpsec = @(v) v * 1.495978e8 / sqrt(1.495978e8^3 / 1.327e11); % [km/s]
kmpsec2emos = @(v) v * sqrt(1.495978e8^3 / 1.327e11) / 1.495978e8; % [EMOS]
au2km = @(d) d * 1.495978e8; % [km]
km2au = @(d) d / 1.495978e8; % [km]

deltav_1_fn = @(mu,r_1,R) ...
    sqrt(mu / r_1) * (sqrt(2 - 2 / (1 + R)) - 1); % [velocity]
deltav_2_fn = @(mu,r_2,R) ...
    sqrt(mu / r_2) * (1 - sqrt(2 - 2 * R / (1 + R))); % [velocity]

%% Given
v_Earth = emos2kmpsec(1); % [km/s]
mu_Sun = 1.327e11; % [km^3/s^2]
mu_Earth = 3.986e5; % [km^3/s^2]
r_Earth = au2km(1); % [km]
r_Venus = au2km(0.7233); % [km]
parking_alt_Earth = 200; % [km]
Earth_radius = 6.37812e3; % [km]
r_park = Earth_radius + parking_alt_Earth; % [km]

%% Calculations
v_a = ...
    sqrt(mu_Sun * (2 / r_Earth - 1 / ((r_Earth + r_Venus) / 2))); % [km/s]

v_inf = v_a - v_Earth; % [km/s]

v_p_hyper = sqrt(v_inf^2 + 2 * mu_Earth / r_park); % [km/s]

v_park = sqrt(mu_Earth / r_park); % [km/s]

deltav = v_p_hyper - v_park; % [km/s]

%% Output
fprintf( ...
    "deltav_1 = %g [km/s]\n", ...
    deltav);