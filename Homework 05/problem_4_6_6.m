% Spring 2024 AER E 351 Homework 05 Problem 4.)6.6)
% Matthew Mehrtens
clear; clc; close all;

%% Given
m_0 = 15000; % [units]
m_L = 1000; % [units]
c = 3048; % [m/s]
m_s = 2000; % [units]

%% Part a.)
m_02 = (m_0^2 * m_L)^(1 / 3); % [units]
m_03 = m_02^2 / m_0; % [units]
lambda = m_02 / (m_0 - m_02); % [units]

%% Part b.)
syms m_s1 m_s2 m_s3;
eqn1 = m_03 * m_s1 - m_L * m_s1 == m_0 * m_s3 - m_02 * m_s3;
eqn2 = m_02 * m_s1 - m_03 * m_s1 == m_0 * m_s2 - m_02 * m_s2;
eqn3 = m_s1 + m_s2 + m_s3 == m_s;
sol = solve([eqn1,eqn2,eqn3],[m_s1,m_s2,m_s3]);
m_s1 = double(sol.m_s1); % [units]
m_s2 = double(sol.m_s2); % [units]
m_s3 = double(sol.m_s3); % [units]

epsilon = m_s1 / (m_0 - m_02); % []

%% Part c.)
syms m_p1 m_p2 m_p3;
eqn1 = m_0 == m_s1 + m_p1 + m_s2 + m_p2 + m_s3 + m_p3 + m_L;
eqn2 = m_02 == m_s2 + m_p2 + m_s3 + m_p3 + m_L;
eqn3 = m_03 == m_s3 + m_p3 + m_L;
sol = solve([eqn1,eqn2,eqn3],[m_p1,m_p2,m_p3]);
m_p1 = double(sol.m_p1); % [units]
m_p2 = double(sol.m_p2); % [units]
m_p3 = double(sol.m_p3); % [units]

%% Part d.)
Z = (1 + lambda) / (epsilon + lambda); % []
deltav = c * log(Z); % [m/s]

%% Output
fprintf( ...
    "Problem 4.)6.6.)a.)\n" + ...
    "m_02 = %g units\n" + ...
    "m_03 = %g units\n" + ...
    "lambda = %g\n\n", ...
    m_02, m_03, lambda);
fprintf( ...
    "Problem 4.)6.6.)b.)\n" + ...
    "m_s1 = %g units\n" + ...
    "m_s2 = %g units\n" + ...
    "m_s3 = %g units\n" + ...
    "epsilon = %g\n\n", ...
    m_s1, m_s2, m_s3, epsilon);
fprintf( ...
    "Problem 4.)6.6.)c.)\n" + ...
    "m_p1 = %g units\n" + ...
    "m_p2 = %g units\n" + ...
    "m_p3 = %g units\n\n", ...
    m_p1, m_p2, m_p3);
fprintf( ...
    "Problem 4.)6.6.)d.)\n" + ...
    "Z = %g\n" + ...
    "deltav = %g m/s\n", ...
    Z, deltav);