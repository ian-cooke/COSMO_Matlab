%MS2
% Milestone 2 Script
% Numerical Integrator

% RK4 Numerical Integrator
% State: X = [sigma, omega]
% This function can be found at the very bottom of this file.

% Consistent tstep for integration and moment of inertia tensor
global tstep;
global I;
global gain;

% Integrate the state with my integrator for 100 seconds
sigma0 = [0.3, 0.2, 0.4]';
omega0 = deg2rad([15, 8, 12]'); % [rad/s]
x0 = [sigma0; omega0];
u.name = 'none';
tspan = 0:tstep:100; % [s]

[state, ~, ~] = intmrprk4(tspan, I, x0, u, 0);
sigma = state(1:3, :);
omega = state(4:6, :);

f1 = figure(1);
set(f1, 'defaultaxesfontsize', 16)
set(f1, 'Visible', 'on')
subplot(2,1,1)
xlabel('Time t [s]')
ylabel('MRP \sigma')
hold on
grid minor
plot(tspan, sigma(1, :), 'Linewidth', 3)
plot(tspan, sigma(2, :), 'Linewidth', 3)
plot(tspan, sigma(3, :), 'Linewidth', 3)
leg = legend('\sigma_1', '\sigma_2', '\sigma_3');
set(leg, 'Location', 'best')
hold off

subplot(2,1,2)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega')
hold on
grid minor
plot(tspan, rad2deg(omega(1, :)), 'Linewidth', 3)
plot(tspan, rad2deg(omega(2, :)), 'Linewidth', 3)
plot(tspan, rad2deg(omega(3, :)), 'Linewidth', 3)
leg = legend('\omega_1', '\omega_2', '\omega_3');
set(leg, 'Location', 'best')
hold off
% Plot Angular Momentum and Kinetic Energy (should be conserved)

H = zeros(1, length(tspan));
KE = H;
for i = 1:length(tspan)
    Iomeg = I*omega(:, i);
    H(i) = norm(Iomeg);
    KE(i) = 1/2*norm(dot(Iomeg, omega(:, i)));
end

f3 = figure(3);
set(f3, 'defaultaxesfontsize', 16)
set(f3, 'Visible', 'on')
xlabel('Time t [s]')
ylabel('Angular Momentum H [kg-m^2/s]')
hold on
yyaxis left
plot(tspan, H, 'Linewidth', 3)
ylim([2.0 2.1])
yyaxis right
plot(tspan, KE, 'Linewidth', 3)
ylabel('Kinetic Energy T [J]')
ylim([0.3 0.4])
hold off
% The angular momentum magnitude and kinetic energy are indeed constant.