%MS5
% Milestone 5 Script
% Orbit Inclinations v/s control performance

% Global vars
global inc;
global tstep;
global T;

% Case 1 - 15 deg inclination and modulating control
inc = deg2rad(15);
% Re-calculate the magnetic field for new inclination
tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
Hb = calcMagField(tspan);

% Now integrate
sigma0 = [0.3,0.2,0.4]'; % [none]
omega0 = deg2rad([15, 8, 12]'); % [rad/s]
x0 = [sigma0;omega0];
control = 'mod';
[Xcase7, uCase7, mCase7] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase7 = Xcase7(1:3, :);
omegaCase7 = Xcase7(4:6, :);

% Only need to plot control error
f18 = figure(18);
set(f18, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [rad/s]')
hold on
grid minor
plot(tspan, rad2deg(omegaCase7(1, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase7(2, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase7(3, :)), 'Linewidth', 2)
% threshold
plot(tspan, 3*ones(1,length(tspan)),'k')
plot(tspan, -3*ones(1,length(tspan)),'k')
leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
set(leg, 'Location', 'best')

% Case 2 - 15 deg inclination and bang-bang control
control = 'bang';
[Xcase8, uCase8, mCase8] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase8 = Xcase8(1:3, :);
omegaCase8 = Xcase8(4:6, :);

% Only need to plot control error
f19 = figure(19);
set(f19, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [rad/s]')
hold on
grid minor
plot(tspan, rad2deg(omegaCase8(1, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase8(2, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase8(3, :)), 'Linewidth', 2)
% threshold
plot(tspan, 3*ones(1,length(tspan)),'k')
plot(tspan, -3*ones(1,length(tspan)),'k')
leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
set(leg, 'Location', 'best')

% Figure out minimum dipole needed for each control for 15 deg inclination
inc = deg2rad(15);
thresh = deg2rad(1); % 1 deg/s threshold


% Case 3 - 105 deg inclination with modulating control
inc = deg2rad(105);
% Re-calculate the magnetic field for new inclination
tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
Hb = calcMagField(tspan);

% Now integrate
control = 'mod';
[Xcase9, uCase9, mCase9] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase9 = Xcase9(1:3, :);
omegaCase9 = Xcase9(4:6, :);

% Only need to plot control error
f20 = figure(20);
set(f20, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [rad/s]')
hold on
grid minor
plot(tspan, rad2deg(omegaCase9(1, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase9(2, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase9(3, :)), 'Linewidth', 2)
% threshold
plot(tspan, 3*ones(1,length(tspan)),'k')
plot(tspan, -3*ones(1,length(tspan)),'k')
leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
set(leg, 'Location', 'best')

% Case 4 - 105 deg inclination with bang-bang control
control = 'bang';
[Xcase10, uCase10, mCase10] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase10 = Xcase10(1:3, :);
omegaCase10 = Xcase10(4:6, :);

% Only need to plot control error
f21 = figure(21);
set(f21, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [rad/s]')
hold on
grid minor
plot(tspan, rad2deg(omegaCase10(1, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase10(2, :)), 'Linewidth', 2)
plot(tspan, rad2deg(omegaCase10(3, :)), 'Linewidth', 2)
% threshold
plot(tspan, 3*ones(1,length(tspan)),'k')
plot(tspan, -3*ones(1,length(tspan)),'k')
leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
set(leg, 'Location', 'best')