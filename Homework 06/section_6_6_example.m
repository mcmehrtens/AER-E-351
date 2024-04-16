% Spring 2024 AER E 351 Homework 06 Problem 1.)6.11.)
% Matthew Mehrtens
clear; clc; close all;

%% Given
m_L = 1000; % [units]
deltaV_T = 10.42; % [km/s]
c = [3.7,3.7,4.2]; % [km/s]
epsilon = [0.08,0.09,0.10]; % []

options = optimset('Display','off');

%% Equations
f = @(epsilon, Z) sum(log(1 - epsilon) + log(Z) - log(1 - epsilon .* Z)); % []
g = @(deltaV_T, c, Z) deltaV_T - sum(c .* log(Z)); % [km/s]

%% Calculations
Z = fmincon(@(X) f(epsilon, X),ones(1,3),[],[],[],[],[],[], ...
    @(X) rocket_constraint(deltaV_T, c, X, g), options); % []

m_0 = m_L * exp(f(epsilon, Z)); % [units]

%% Output
fprintf( ...
    "Section 6.6 Example:\n" + ...
    "m_0 = %g units\n", ...
    m_0);