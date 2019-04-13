%COSMO
% script to do some controls analysis and monte carlos for the COSMO
% cubesat

clear; close all;

globalvars;

global tstep;
global I;
global T;

% Initial position and angular velocity
sigma0 = [0.2, 0.3, 0.4]';
omega0 = deg2rad([5, 5, 5]'); % [rad/s]
x0 = [sigma0; omega0];
u.name = 'loveraTrack1';
u.sigma_RN = [0.0, 0.0, 0.0]';
u.K_sigma = 0.1;
u.K_omega = 0.07;
u.epsilon = 0.1;
tspan = 0:tstep:30*T; % [s]

% Calculate magnetic field
Hb = calcMagField(tspan);

[state, u, m] = intmrprk4(tspan, I, x0, u, Hb);
sigma = state(1:3, :);
omega = state(4:6, :);

f1 = figure(1);
plotState(f1, tspan, sigma, omega);
f2 = figure(2);
plotControl(f2, tspan, u, m);