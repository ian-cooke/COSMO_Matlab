% Solar Beta Angle Less than 15 deg script
clear; close all;
restoredefaultpath;
addpath('/Users/iancooke/Dropbox/CUBoulder/Grad Year 1/Spring/StatOD/Homeworks/HW7/functions')
addpath('/Users/iancooke/Dropbox/CUBoulder/Grad Year 1/Spring/StatOD/Homeworks/HW7/macros')

%% parameters for sim

% orbit params
orbit_radius = 415.28e3;
R_E = 6378e3;
a = (R_E + orbit_radius);
e = 0.0003449;
inc = 51.7;
Omega = 0;
omega = 0;
nu_0 = 0;
inc_sun = 17;
inc_eff = inc - inc_sun;
au = 149597870700;

% convert
mu = 3.986004415e14; % gravitational parameter
[R_0, V_0] = kep2eci(a, e, inc_eff, Omega, omega, nu_0, mu);
x_0 = [R_0; V_0];

% other stuff
n = 6;
dt = 1;
T = 2*pi*a^(3/2)/sqrt(mu);
sun_pos = 0:pi/50:2*pi; % radial sun position
tspan = 0:dt:T;
J2 = 1.082626925638815e-3;
power_peak = 20.65; % [W]

% params struct for ode fxn
params.R_E = R_E;
params.n = n;
params.mu = mu;
params.J2 = J2;


% ode options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 

% integrate the trajectory from 0 to T
[~, x_temp] = ode113(@(t,x) keplerJ2_wPhi_ODE(t, x, params), tspan, [x_0; reshape(eye(n), n^2, 1)], options);
orbit = x_temp(:, 1:n)';

% loop vars
frac_less_than_15 = zeros(length(sun_pos),1);
power_frac = zeros(length(sun_pos), length(tspan));


% loop over sun angle and orbit
for ii = 1:length(sun_pos)
    sun_vec = au.*[cos(sun_pos(ii)), sin(sun_pos(ii)), 0]';
    % loop over orbit
    sum = 0;
    for jj = 1:length(tspan)
        r = orbit(1:3, jj);
        sun_to_sat = sun_vec + r;
        sun_to_sat = sun_to_sat./norm(sun_to_sat);
        r_unit = r./norm(r);
        beta = acos(dot(r_unit, sun_to_sat));
        power_frac(ii, jj) = dot(r_unit, sun_to_sat);
        if power_frac(ii, jj) < 0
            power_frac(ii, jj) = 0;
        end
        if rad2deg(beta) < 15
            sum = sum + 1;
        end
    end
    frac_less_than_15(ii) = sum/length(tspan);
end

power_avg = power_peak.*mean(power_frac, 2);

% plot percentage vs. sun angle
f1 = figure(1);
set(f1, 'defaultaxesfontsize', 16)
set(f1, 'Position', get(0, 'Screensize'));
xlabel('Sun Angle [deg]')
ylabel('Percent [%]')
title('Percentage Of Time In 1 Orbit Solar Panels Are Within 15 deg of Sun Vector vs. Solar Orientation')
hold on
grid on
plot(rad2deg(sun_pos), frac_less_than_15.*100, 'Linewidth', 1.5)

% plot power contour
f2 = figure(2);
set(f2, 'defaultaxesfontsize', 16)
set(f2, 'Position', get(0, 'Screensize'));
xlabel('Sun Orientation [deg]')
ylabel('Spacecraft Mean Anomaly [deg]')
zlabel('Peak Power Fraction, Cosine of Solar Anlge Relative to SC')
title('Peak Power Fraction vs. Solar Orientation and Spacecraft Mean Anomaly')
hold on
grid on
xlim([0 360])
ylim([0 360])
zlim([0 1])
[X, Y] = meshgrid(rad2deg(sun_pos), tspan./tspan(end)*360);
surf(X, Y, power_frac')
shading interp
colorbar

% plot average power vs. sun orientation
f3 = figure(3);
set(f3, 'defaultaxesfontsize', 16)
set(f3, 'Position', get(0, 'Screensize'));
xlabel('Sun Orientation [deg]')
ylabel('Average Spacecraft Power [W]')
title('Average Spacecraft Power Over 1 Orbit vs. Solar Orientation')
hold on
grid on
plot(rad2deg(sun_pos), power_avg, 'Linewidth', 1.5)

% worst case scenario orbit
f4 = figure(4);
set(f4, 'defaultaxesfontsize', 16)
set(f4, 'Position', get(0, 'Screensize'));
xlabel('Time [min]')
ylabel('Spacecraft Power [W]')
title('Worst Case Solar Power Orbit')
hold on
grid on
plot(tspan./60, power_peak.*power_frac(25,:), 'Linewidth', 1.5)
            