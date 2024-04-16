% Spring 2024 AER E 351 Homework 06 Problem 2.)7.2.)
% Matthew Mehrtens
clear; clc; close all;

%% Given
mu_Sun = 1; % [au^3/ctu^2]
r_Earth = 1; % [au]
r_Mars = 1.5237; % [au]
R = r_Mars / r_Earth; % []

cdu = 1.495978e8; % [km]
mu_sun = 1.327e11; % [km^3/s^2]
ctu = sqrt(cdu^3 / mu_sun); % [s]

c = 0.1; % [EMOS]

m_0 = 100; % [units]
epsilon = 1 / 7; % []

%% Equations
deltaV_1_fn = @(mu,r_1,R) ...
    sqrt(mu / r_1) * (sqrt(2 - 2 / (1 + R)) - 1); % [velocity]
deltaV_2_fn = @(mu,r_2,R) ...
    sqrt(mu / r_2) * (1 - sqrt(2 - 2 * R / (1 + R))); % [velocity]

ctu2days = @(t) t * ctu / 86400; % [days]

emos2kmpsec = @(emos) emos * cdu / ctu; % [km/s]

t_H_fn = @(r_1,r_2,mu) pi * sqrt((r_1 + r_2)^3 / (8 * mu)); % [time]

m_pbym_0_fn = @(deltaV,c) 1 - exp(-deltaV / c); % []

m_L_eqn = @(m_0, Z, epsilon) ...
    m_0 * (1 + (Z - 1) / (1 - Z * epsilon))^(-1); % [mass]

%% Calculations
deltaV_1 = deltaV_1_fn(mu_Sun,r_Earth,R); % [EMOS]
deltaV_2 = deltaV_2_fn(mu_Sun,r_Mars,R); % [EMOS]
deltaV = deltaV_1 + deltaV_2; % [EMOS]

t_H = ctu2days(t_H_fn(r_Earth,r_Mars,mu_Sun)); % [days]

m_pbym_0 = m_pbym_0_fn(deltaV,c); % []

m_L = m_L_eqn(m_0,exp(deltaV / c),epsilon); % [units]

%% Output
fprintf( ...
    "Problem 2.)7.2.)\n" + ...
    "a.)\n" + ...
    "deltaV_1 = %g EMOS = %g km/s\n" + ...
    "deltaV_2 = %g EMOS = %g km/s\n" + ...
    "deltaV = %g EMOS = %g km/s\n\n", ...
    deltaV_1,emos2kmpsec(deltaV_1),deltaV_2,emos2kmpsec(deltaV_2), ...
    deltaV,emos2kmpsec(deltaV));
fprintf( ...
    "b.)\n" + ...
    "t_H = %g days\n\n", ...
    t_H)
fprintf( ...
    "c.)\n" + ...
    "m_p / m_0 = %g\n\n", ...
    m_pbym_0);
fprintf( ...
    "d.)\n" + ...
    "m_L = %g units\n", ...
    m_L);