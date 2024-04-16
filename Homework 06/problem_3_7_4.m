% Spring 2024 AER E 351 Homework 06 Problem 3.)7.4.)
% Matthew Mehrtens
clear; clc; close all;

%% Given
r_Earth = 6.37812e3; % [km]
r_1 = r_Earth + 278; % [km]

mu_Earth = 3.986e5; % [km^3/s^2]

T_2 = 23.934; % [hours]
T_2 = T_2 * 3600; % [s]

theta = deg2rad(28); % [rad]

g = 9.80665; % [m/s^2]

I_sp = 300; % [s]
m_0 = 100; % [units]
epsilon = 1 / 7; % []

%% Equations
r_c_fn = @(T,mu) (T^2 * mu / (4 * pi^2))^(1 / 3);

t_H_fn = @(r_1,r_2,mu) pi * sqrt((r_1 + r_2)^3 / (8 * mu)); % [time]

v_c_fn = @(mu,r) sqrt(mu / r); % [velocity]
v_p_fn = @(mu,r_1,R) ...
    sqrt(2 * mu / r_1) * sqrt(1 - 1 / (1 + R)); % [velocity]
v_A_fn = @(mu,r_2,R) ...
    sqrt(2 * mu / r_2) * sqrt(1 - R / (1 + R)); % [velocity]

deltav_1_fn = @(v_p,v_c1) v_p - v_c1; % [velocity]
deltav_2_fn = @(v_A,v_c2,theta) ...
    sqrt(v_A^2 + v_c2^2 - 2 * v_A * v_c2 * cos(theta)); % [velocity]

m_L_fn = @(m_0,Z,epsilon) ...
    m_0 * (1 + (Z - 1) / (1 - Z * epsilon))^(-1); % [mass]
m_s_fn = @(m_0,Z,m_L) m_0 / Z - m_L; % [mass]
m_p_fn = @(m_0,m_s,m_L) m_0 - m_s - m_L; % [mass]

%% Calculations
% a.)
r_2 = r_c_fn(T_2,mu_Earth); % [km]

t_H = t_H_fn(r_1,r_2,mu_Earth); % [s]
t_H = t_H / 3600; % [hours]

% b.)
R = r_2 / r_1; % []
v_c1 = v_c_fn(mu_Earth,r_1); % [km/s]
v_c2 = v_c_fn(mu_Earth,r_2); % [km/s]
v_p = v_p_fn(mu_Earth,r_1,R); % [km/s]
v_A = v_A_fn(mu_Earth,r_2,R); % [km/s]

deltav_1 = deltav_1_fn(v_p,v_c1); % [km/s]
deltav_2 = deltav_2_fn(v_A,v_c2,theta); % [km/s]

% c.)
c = I_sp * g; % [m/s]
c = c / 1000; % [km/s]

deltav = deltav_1 + deltav_2; % [km/s]

Z = exp(deltav / c); % []
m_L = m_L_fn(m_0,Z,epsilon); % [units]
m_s = m_s_fn(m_0,Z,m_L); % [units]
m_p = m_p_fn(m_0,m_s,m_L); % [units]

%% Output
fprintf( ...
    "Problem 3.)7.4.)\n" + ...
    "a.)\n" + ...
    "t_H = %g hours\n\n", ...
    t_H);
fprintf( ...
    "b.)\n" + ...
    "deltav_1 = %g km/s\n" + ...
    "deltav_2 = %g km/s\n\n", ...
    deltav_1, deltav_2);
fprintf( ...
    "c.)\n" + ...
    "Z = %g\n" + ...
    "m_p = %g units\n" + ...
    "m_s = %g units\n" + ...
    "m_L = %g units\n", ...
    Z, m_p, m_s, m_L);