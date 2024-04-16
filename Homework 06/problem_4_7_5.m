% Spring 2024 AER E 351 Homework 06 Problem 4.)7.5.)
% Matthew Mehrtens
clear; clc; close all;

%% Given
r_1 = 1; % [au]
r_2 = 1.524; % [au]
R = r_2 / r_1; % []

T_Earth = 365.256; % [days]
T_Earth = T_Earth * 86400; % [s]

T_Mars = 365.256 + 321.73; % [days]
T_Mars = T_Mars * 86400; % [s]

mu = 1; % [au^3/ctu^2]

au = 1.495978e8; % [km]
mu_Sun = 1.327e11; % [km^3/s^2]
ctu = sqrt(au^3 / mu_Sun); % [s]

S_Mars = 778; % [days]
S_Mars = S_Mars * 86400; % [s]

%% Equations
beta_12_fn = @(R) pi * (1 - ((1 + R) / (2 * R))^(3 / 2)); % [rad]
beta_21_fn = @(R) pi * (1 - ((1 + R) / 2)^(3 / 2)); % [rad]

n_fn = @(T) 2 * pi / T; % [rad/time]

t_H_fn = @(r_1,r_2,mu) pi * sqrt((r_1 + r_2)^3 / (8 * mu)); % [time]

ctu2sec = @(t) t * ctu; % [s]

t_opp_fn = @(beta_12,n_Earth,n_Mars) beta_12 / (n_Earth - n_Mars); % [time]

t_transfer_2_fn = @(beta_21,theta,theta_2,n_Earth,n_Mars) ...
    (beta_21 + theta - theta_2) / (n_Earth - n_Mars); % [time]

%% Calculations
% a.)
beta_12 = beta_12_fn(R); % [rad]
beta_21 = beta_21_fn(R); % [rad]

% b.)
n_Earth = n_fn(T_Earth); % [rad/s]
t_H = t_H_fn(r_1,r_2,mu); % [ctu]
t_H = ctu2sec(t_H); % [s]
theta_2 = n_Earth * t_H; % [rad]255

% c.)
n_Mars = n_fn(T_Mars); % [rad/s]
t_opp = t_opp_fn(beta_12,n_Earth,n_Mars); % [s]
t_opp = t_opp / 86400; % [days]

% d.)
theta_opp = mod(n_Earth * S_Mars,2 * pi); % [rad]

% e.)
t_transfer_2 = t_transfer_2_fn(beta_21,pi,theta_2,n_Earth,n_Mars) ...
    + S_Mars; % [s]
t_transfer_2 = t_transfer_2 / 86400; % [days]

%% Output
fprintf( ...
    "Problem 4.)7.5.)\n" + ...
    "a.)\n" + ...
    "beta_12 = %g rad = %g°\n" + ...
    "beta_21 = %g rad = %g°\n\n", ...
    beta_12,rad2deg(beta_12),beta_21,rad2deg(beta_21));
fprintf( ...
    "b.)\n" + ...
    "n_Earth = %g rad/s\n" + ...
    "t_H = %g s\n" + ...
    "theta_2 = %g rad = %g°\n\n", ...
    n_Earth,t_H,theta_2,rad2deg(theta_2));
fprintf( ...
    "c.)\n" + ...
    "n_Mars = %g rad/s\n" + ...
    "t_opp = %g days\n\n", ...
    n_Mars, t_opp);
fprintf( ...
    "d.)\n" + ...
    "theta_opp = %g rad = %g°\n\n", ...
    theta_opp, rad2deg(theta_opp));
fprintf( ...
    "e.)\n" + ...
    "t_transfer_2 = %g days\n", ...
    t_transfer_2);
