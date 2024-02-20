% AER E 351 Homework 02 Problem 2.6 Part d
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Given
R_Earth = 6.37812e3; % [km]
mu_Earth = 3.986e5; % [km^3/s^2]

%% Calculations
R = 144 * R_Earth; % [km]
delta_t = sqrt(2 / 9 * R / mu_Earth) * R; % [s]
delta_t_days = delta_t / 3600 / 24; % [days]
total_time = delta_t_days * 2; % [days]

fprintf( ...
    "delta_t = %g s\n" + ...
    "delta_t = %g days\n" + ...
    "total_time = %g days\n", ...
    delta_t, delta_t_days, total_time);