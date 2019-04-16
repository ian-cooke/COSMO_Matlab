%COSMO
% script to do some controls analysis and monte carlos for the COSMO
% cubesat

clear; close all;

globalvars;

global tstep;
global I;
global T;

% Initial state & other params
sigma0 = [0.2, 0.3, 0.4]';
omega0 = deg2rad([0, 0, 0]'); % [rad/s]
x0 = [sigma0; omega0];
u.name = 'bdotTrack';
u.sigma_RN = [0.5, 0.0, -0.5]';
u.K_sigma = 0.00000001;
u.K_omega = 0.5;
u.epsilon = 0.005;
tspan = 0:tstep:15*T; % [s]

% Calculate magnetic field
Hb = calcMagField(tspan);

[state, u, m] = intmrprk4(tspan, I, x0, u, Hb);
sigma = state(1:3, :);
omega = state(4:6, :);

f1 = figure(1);
plotState(f1, tspan, sigma, omega);
f2 = figure(2);
plotControl(f2, tspan, u, m);