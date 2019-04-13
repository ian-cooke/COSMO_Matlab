%MINM
% Figure out minimum dipole moments for each control law to get omega below
% 1 deg/s threshold at end of 3 orbits for 15 deg inclination case
% Global vars
global inc;
global tstep;
global T;
global m_max;

% Case 1 - mod control, init conditions mission overview
inc = deg2rad(15);
% Re-calculate the magnetic field for new inclination
tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
Hb = calcMagField(tspan);
m_max = 4.33;
% Now integrate
sigma0 = [0.3,0.2,0.4]'; % [none]
omega0 = deg2rad([15, 8, 12]'); % [rad/s]
x0 = [sigma0;omega0];
% control = 'mod';
% [X, u, m] = intmrprk4(tspan, I, x0, control, Hb);
% sigma = X(1:3, :);
% omega = X(4:6, :);
% 
% % Only need to plot control error
% f18 = figure(18);
% set(f18, 'defaultaxesfontsize', 16)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega [rad/s]')
% hold on
% grid minor
% plot(tspan, rad2deg(omega(1, :)), 'Linewidth', 2)
% plot(tspan, rad2deg(omega(2, :)), 'Linewidth', 2)
% plot(tspan, rad2deg(omega(3, :)), 'Linewidth', 2)
% % threshold
% plot(tspan, 1*ones(1,length(tspan)),'k')
% plot(tspan, -1*ones(1,length(tspan)),'k')
% leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
% set(leg, 'Location', 'best')

% Case 2 - bang control, init conditions mission overview
% Re-calculate the magnetic field for new inclination
m_max = 4.15;
control = 'bang';
[X, u, m] = intmrprk4(tspan, I, x0, control, Hb);
sigma = X(1:3, :);
omega = X(4:6, :);

% Only need to plot control error
f18 = figure(18);
set(f18, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [rad/s]')
hold on
grid minor
plot(tspan, rad2deg(omega(1, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omega(2, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omega(3, :)), 'Linewidth', 2)
% threshold
plot(tspan, 1*ones(1,length(tspan)),'k')
plot(tspan, -1*ones(1,length(tspan)),'k')
leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
set(leg, 'Location', 'best')