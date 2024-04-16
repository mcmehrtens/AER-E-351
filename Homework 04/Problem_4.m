% Spring 2024 AER E 351 Homework 04 Problem 4
% Matthew Mehrtens
clear, clc, close all;

%% Given
r_1 = 1; % [au]
r_2 = 1.524; % [au]
theta = deg2rad(75); % [rad]

t_F_a = 0.17; % [years]
t_F_b = 0.34; % [years]
t_F_c = 0.90; % [years]

mu = 1; % [au^3/ctu^2]

cdu = 1.495978e8; % [km]
mu_sun = 1.327e11; % [km^3/s^2]
ctu = sqrt(cdu^3 / mu_sun); % [s]

Earth_siderial = 365.256; % [days]

ctu2years = @(t) t * ctu / 86400 / Earth_siderial; % [years]
years2ctu = @(t) t * 86400 * Earth_siderial / ctu; % [ctu]

%% Equations
c_fn = @(r_1, r_2, theta) sqrt(r_1^2 + r_2^2 ...
    - 2 * r_1 * r_2 * cos(theta)); % [distance]
s_fn = @(r_1, r_2, c) (r_1 + r_2 + c) / 2; % [distance]
t_p_fn = @(mu, s, c, theta) sqrt(2) / (3 * mu) * (s^(3 / 2) ...
    - sign(sin(theta)) * (s - c)^(3 / 2)); % [time]
lamberts_eqn = @(mu, t_F, s, c, a) -sqrt(mu) * t_F ...
    + a^(3 / 2) * (2 * asin(s / (2 * a)) - 2 * asin((s - c) / (2 * a)) ...
    - (sin(2 * asin(sqrt(s / (2 * a)))) ...
    - sin(2 * asin(sqrt((s - c) / (2 * a))))));
lamberts_eqnsharp = @(mu, t_F, s, c, a) -sqrt(mu) * t_F ...
    + a^(3 / 2) ...
    * (2 * pi - 2 * asin(s / (2 * a)) - 2 * asin((s - c) / (2 * a)) ...
    - (sin(2 * pi - 2 * asin(sqrt(s / (2 * a)))) ...
    - sin(2 * asin(sqrt((s - c) / (2 * a))))));

%% Common Calculations
c = c_fn(r_1, r_2, theta); % [au]
s = s_fn(r_1, r_2, c); % [au]
t_p = t_p_fn(mu, s, c, theta); % [ctu]

alpha_m = pi; % [rad]
beta_m = 2 * asin(sqrt((s - c) / s)); % [rad]
t_m = sqrt(s^3 / (8 * mu)) * (pi - beta_m + sin(beta_m)); % [ctu]

%% Part b.)          
a_b = fsolve(@(a) lamberts_eqn(mu, years2ctu(t_F_b), s, c, a), ...
    1); % [au]

%% Part c.)
a_c = fsolve(@(a) lamberts_eqnsharp(mu, years2ctu(t_F_c), s, c, a), ...
    1); % [au]

%% Output
fprintf(...
    "Problem 4.)\n" + ...
    "Given:\n" + ...
    "r_1 = %g au\n" + ...
    "r_2 = %g au\n" + ...
    "theta = %g° = %g rad\n" + ...
    "mu = %g au^3/ctu^2\n" + ...
    "t_F_a = %g years\n" + ...
    "t_F_b = %g years\n" + ...
    "t_F_c = %g years\n\n" + ...
    "Common Calculations:\n" + ...
    "c = %g au\n" + ...
    "s = %g au\n" + ...
    "t_p = %g ctu = %g years\n" + ...
    "t_m = %g ctu = %g years\n\n", ...
    r_1, r_2, rad2deg(theta), theta, mu, t_F_a, t_F_b, t_F_c, c, s, ...
    t_p, ctu2years(t_p), t_m, ctu2years(t_m));
fprintf(...
    "Problem 4.)a.)\n" + ...
    "t_F_a < t_p; so, no elliptic transfer orbit exists.\n\n");
fprintf(...
    "Problem 4.)b.)\n" + ...
    "a_b = %g au\n\n", ...
    a_b);
fprintf(...
    "Problem 4.)c.)\n" + ...
    "a_c = %g au\n", ...
    a_c);