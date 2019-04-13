%MS3
% Milestone 3 Script
% Geomagnetic Modeling
%% Milestone 3: Geomagnetic Modeling
% Modeling the tilted magnetic dipole using the coordinate frame representations

% Orbital parameters
global h;

% LVLH reference frame angles
global Omega;
global theta0;
global thetadot;
global T;
global inc;

global tstep;

% Plot the magnetic field  as seen by the orbit frame H over one orbit
% No rk4 needed - no differential equation - simply nee a loop
tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
% Propagate the state
Hb = calcMagField(tspan);

f4 = figure(4);
set(f4, 'Visible', 'on')
set(f4, 'defaultaxesfontsize', 16)
hold on
grid minor
xlabel('Time t [s]')
xlim([0 T])
ylabel('Magnetic Field Strength [T]')
plot(tspan, Hb(1, :), 'Linewidth', 3)
plot(tspan, Hb(2, :), 'Linewidth', 3)
plot(tspan, Hb(3, :), 'Linewidth', 3)
leg = legend('^Hb_1', '^Hb_2', '^Hb_3');
set(leg, 'Location', 'best')

%Plot the Magnetic Field as seen by the body frame for 100 seconds
% Just need to convert H coordinates to B coordinates
% Bb = [BN][HN]'*Hb
tspan = 0:tstep:100;
h = tspan(2) - tspan(1);
t = tspan(1);
ind = 1;
while (t <= tspan(end))
    % [BN] DCM using MRPs calculated in Milestone 2
    s = sigma(:, ind);
    sq = s'*s;
    BN = eye(3) + (8*tilde(s)*tilde(s) - 4*(1-sq)*tilde(s))/(1+sq)^2;
                 
    % [HN] DCM
    theta = theta0 + thetadot*t;
    HN = [cos(Omega)*cos(theta) - sin(Omega)*cos(inc)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(inc)*sin(theta), sin(inc)*sin(theta);
          -cos(Omega)*sin(theta) - sin(Omega)*cos(inc)*cos(theta), cos(Omega)*cos(inc)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(inc);
          sin(Omega)*sin(inc), -cos(Omega)*sin(inc), cos(inc)];
      
    % [BH] DCM
    BH = BN*HN';
    Bb(:, ind) = BH*Hb(:, ind);
    
t = t + h;
ind = ind + 1;
end

% Plot it
f5 = figure(5);
set(f5, 'Visible', 'on')
set(f5, 'defaultaxesfontsize', 16)
hold on
grid minor
xlabel('Time t [s]')
ylabel('Magnetic Field Strength [T]')
plot(tspan, Bb(1, :), 'Linewidth', 3)
plot(tspan, Bb(2, :), 'Linewidth', 3)
plot(tspan, Bb(3, :), 'Linewidth', 3)
leg = legend('^Bb_1', '^Bb_2', '^Bb_3');
set(leg, 'Location', 'best')
