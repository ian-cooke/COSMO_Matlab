%COSMO_fullstate
% script to do some controls analysis and monte carlos for the COSMO
% cubesat

clear; close all;

globalvars;

global tstep;
global I;
global T;
global r;
global mu;
global inc;
global ecc;
global theta0;
global Omega;
global arg_per;

% Initial state & other params
sigma0 = [0.2, 0.3, 0.4]';
omega0 = deg2rad([1, 1, 1]'); % [rad/s]
% orbit elements
[pos0, vel0] = kep2eci(r, ecc, rad2deg(inc), Omega, arg_per, rad2deg(theta0), mu);
x0 = [sigma0; omega0; pos0; vel0];
u.name = 'loveraTrack1';
u.sigma_RN = [0.5, 0.0, -0.5]';
u.K_sigma = 0.0001;
u.K_omega = 50;
u.epsilon = 0.005;
u.use_velocitypoint = 0;
tspan = 0:tstep:15*T; % [s]

% Calculate magnetic field
Hb = calcMagField(tspan);

%[state, u, m] = intmrprk4(tspan, I, x0, u, Hb);
[state, u, m, params] = int_mrp_orbit_rk4(tspan, I, x0, u, Hb);
sigma = state(1:3, :);
omega = state(4:6, :);
pos = state(7:9, :);

f1 = figure(1);
plotState(f1, tspan, sigma, omega);
f2 = figure(2);
plotControl(f2, tspan, u, m);

f3 = figure(3);
hold on
plot(tspan, params.sigma_BR(1, :))
plot(tspan, params.sigma_BR(2, :))
plot(tspan, params.sigma_BR(3, :))
