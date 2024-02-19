% AER E 351 Homework 02 Problem 2.13
% Lucy Gates and Matthew Mehrtens
clear, clc, close all;

%% Given
R_Earth = 6.378e3; % [km]
mu_Earth = 3.986e5; % [km^3/s^2]
au = 1.495978e8; % [km]
mu_sun = 1.327e11; % [km^3/s^2]

%% Calculations b
CTU = sqrt(R_Earth^3 / mu_Earth); % [s]
CTU_minutes = CTU / 60; % [min]

%% Calculations c
CDU = mu_Earth / (R_Earth)^2; % [km]
CDU_meters = CDU * 1000; % [m]

%% Calculations d
CTU_sun = sqrt(au^3 / mu_sun); % [s]
CTU_sun_days = CTU_sun / 3600 / 24; % [days]

fprintf( ...
    "(b) CTU = %g min\n" + ...
    "(c) CDU = %g m\n" + ...
    "(d) CTU = %g days\n", ...
    CTU_minutes, CDU_meters, CTU_sun_days);