% AER E 351 Homework 02 Problem 2.6 Part d
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Given
R_Earth = 6.37812e3; % [km]
mu_Earth = 3.986e5; % [km^3/s^2]

%% Calculations
delta_t = sqrt(288 * R_Earth^3 / mu_Earth); % [s]
delta_t_days = delta_t / 3600 / 24; % [days]

fprintf( ...
    "delta_t = %g s\n" + ...
    "delta_t = %g days\n", ...
    delta_t, delta_t_days);