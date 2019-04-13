% COSMO gravity gradient stabilization analysis
globalvars;

global I;
global T;
global thetadot;

sigma0 = [0.4,0.2,0.3]';
omega0 = deg2rad([0,0,0]');
x0 = [sigma0;omega0];

tspan = 0:1:T*10;

x = intmrprk4_gg(tspan, I, x0);
sigma = x(1:3,:);
omega = x(4:6,:);

f1 = figure(1);
plotState(f1, tspan, sigma, omega)

% do it with EOMs
% state - [s1,s2,s3,s1dot,s2dot,s3dot]'
A = [zeros(3),eye(3);...
     2*thetadot^2*(I3-I2)/I1, 0,0,0,0, thetadot*((I3-I2)/I1-1);...
     0, 2*thetadot^2*(I1-I3)/I2, 0, 0, 0, 0;...
     0,0, thetadot^2*(I2-I1)/I3, 0,0, thetadot*((I2-I1)/I3-1)];

dt = 1;
F = expm(A*dt);
state0 = [0.4,0.3,0.2,deg2rad(0),0,deg2rad(0)];
state = zeros(6, length(tspan));
state(:,1) = state0;
for i = 2:length(tspan)
    state(:, i) = F*state(:, i-1);
    if (norm(state(1:3,i)) > 1)
        state(1:3,i) = -state(1:3,i)/norm(state(1:3,i))^2;
    end
end

f2 = figure(2);
subplot(2,1,1)
hold on
plot(tspan, state(1,:))
plot(tspan, state(2,:))
plot(tspan, state(3,:))
subplot(2,1,2)
hold on
plot(tspan, state(4,:))
plot(tspan, state(5,:))
plot(tspan, state(6,:))
