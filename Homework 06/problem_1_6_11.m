% Spring 2024 AER E 351 Homework 06 Problem 1.)6.11.)
% Matthew Mehrtens
clear; clc; close all;

%% Given
m_L = 1000; % [kg]
deltaV_T = 8.5; % [km/s]
deltaV_T = deltaV_T * 1000; % [m/s]
c = [3500,3800,4100]; % [m/s]
epsilon = [0.10,0.12,0.09]; % []
gravity = 9.81; % [m/s^2]
t = 90; % [s]

options = optimset('Display','off');

%% Equations
f = @(epsilon, Z) ...
    sum(log(1 - epsilon) + log(Z) - log(1 - epsilon .* Z)); % []
g = @(deltaV_T, c, Z) deltaV_T - sum(c .* log(Z)); % [m/s]
g_gravity = @(deltaV_T, c, Z, gravity, t) g(deltaV_T, c, Z) + gravity * t; % [m/s]

%% Calculations
Z = fmincon(@(X) f(epsilon, X),ones(1,3),[],[],[],[],[],[], ...
    @(X) rocket_constraint(deltaV_T, c, X, g), options); % []

m_0 = m_L * exp(f(epsilon, Z)); % [kg]

Z_gravity = fmincon(@(X) f(epsilon, X),ones(1,3),[],[],[],[],[],[], ...
    @(X) rocket_constraint( ...
    deltaV_T, c, X, @(deltaV_T,c,X) g_gravity(deltaV_T,c,X,gravity,t)), ...
    options); % []

m_0_gravity = m_L * exp(f(epsilon, Z_gravity)); % [kg]

%% Output
fprintf( ...
    "Problem 1.)6.11.)\n" + ...
    "a.)\n" + ...
    "m_0 = %g kg\n" + ...
    "b.)\n" + ...
    "m_0 = %g kg\n", ...
    m_0, m_0_gravity);