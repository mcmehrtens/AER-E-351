% Spring 2024 AER E 351 Section 9.3 Example (earth—Venus Hohmann transfer)
% Matthew Mehrtens
clear; close all; clc;

%% Given
G = 6.67259e-11 / 1000^3; % [km^3/(kg*s^2)]
m_Earth = 5.974e24; % [kg]
m_Venus = 0.8149 * m_Earth; % [kg]
mu_Earth = G * m_Earth; % [km^3/s^2]
mu_Venus = G * m_Venus; % [km^3/s^2]
mu_Sun = 1.327e11; % [km^3/s^2]
alt_Earth = 200; % [km]
alt_Venus = 500; % [km]
r_Earth = 6.37812e3; % [km]
r_Venus = 0.949 * r_Earth; % [km]
r_c1 = r_Venus + alt_Venus; % [km]
r_c2 = r_Earth + alt_Earth; % [km]
r_1 = 0.7233; % [au]
r_2 = 1; % [au]

cdu = 1.495978e8; % [km]
ctu = sqrt(cdu^3 / mu_Sun); % [s]

%% Equations
v_c_fn = @(mu,r) sqrt(mu / r); % [velocity]

deltav_1_fn = @(mu,r_1,R) ...
    sqrt(mu / r_1) * (sqrt(2 - 2 / (1 + R)) - 1); % [velocity]
deltav_2_fn = @(mu,r_2,R) ...
    sqrt(mu / r_2) * (1 - sqrt(2 - 2 * R / (1 + R))); % [velocity]

au2km = @(au) au * cdu; % [km]
emos2kmpsec = @(emos) emos * cdu / ctu; % [km/s]
kmpsec2emos = @(kmpsec) kmpsec * ctu / cdu; % [EMOS]

%% Calculations
R = r_2 / r_1; % []
v_c1 = v_c_fn(mu_Venus,r_c1); % [km/s]
v_c2 = v_c_fn(mu_Earth,r_c2); % [km/s]
deltav_1 = deltav_1_fn(mu_Sun,au2km(r_1),R); % [km/s]
deltav_2 = deltav_2_fn(mu_Sun,au2km(r_2),R); % [km/s]

v_pe = sqrt(2 * mu_Earth / r_c2 + deltav_2^2) % [km/s]

%% Output
fprintf( ...
    "R = %g\n" + ...
    "v_c1 = %g km/s\n" + ...
    "v_c2 = %g km/s\n" + ...
    "deltav_1 = %g km/s\n" + ...
    "deltav_2 = %g km/s\n" + ...
    "deltav_2 / v_c1 = %g\n", ...
    R, v_c1, v_c2, deltav_1, deltav_2, deltav_2 / v_c1);