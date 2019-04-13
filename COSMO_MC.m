%COSMO_MC
% Monte Carlo Analysis

% Global variables
globalvars;
global T;
global tstep;
global I;
global inc;

% Modulating control
N = 25; % For 25 iterations
tspan = 0:tstep:3*T;
sigma0 = [0.3, 0.2, 0.4]';
control = 'mod';
inc = deg2rad(45);
Hb = calcMagField(tspan);

% Lower and upper bounds
lower = 4.5;
upper = 5.5;

% Figure set-up
f22 = figure(22);
set(f22, 'defaultaxesfontsize',16)


for ind = 1:N
    % Step 1 - create random angular velocity vector
    omega_rand0 = (upper - lower).*rand(3, 1) + lower;
    omega_rand0 = deg2rad(omega_rand0);
    
    % Step 2 - Simulate the random angular velocity
    x0 = [sigma0;omega_rand0];
    [X, ~, m1] = intmrprk4(tspan, I, x0, control, Hb);
    omega = X(4:6, :);
    
    % Step 3 - Calculate the norm and plot
    norm = vecnorm(omega,2,1); % returns euclidean norm (2) of columns (1)
    semilogy(tspan./T, rad2deg(norm), 'Linewidth', 2);
    
    % Step 4 - repeat this 25 times (happening due to loop)
    if ind == 1
        hold on
        xlabel('Orbit #')
        ylabel('Control Error Norm or Angular Velocity Norm [deg/s]')
        title('Angular Velocity vs. Orbit, Modulating b-dot Law')
        grid minor
        semilogy(tspan./T, 3*ones(1,length(tspan)),'k');
    end
    
end
    

% Bang-Bang control
control = 'bang';

% Figure set-up
f23 = figure(23);
set(f23, 'defaultaxesfontsize',16)


for ind = 1:N
    % Step 1 - create random angular velocity vector
    omega_rand0 = (upper - lower).*rand(3, 1) + lower;
    omega_rand0 = deg2rad(omega_rand0);
    
    % Step 2 - Simulate the random angular velocity
    x0 = [sigma0;omega_rand0];
    [X, ~, m2] = intmrprk4(tspan, I, x0, control, Hb);
    omega = X(4:6, :);
    
    % Step 3 - Calculate the norm and plot
    norm = vecnorm(omega,2,1); % returns euclidean norm (2) of columns (1)
    semilogy(tspan./T, rad2deg(norm), 'Linewidth', 2);
    
    % Step 4 - repeat this 25 times (happening due to loop)
    if ind == 1
        hold on
        xlabel('Orbit #')
        ylabel('Control Error Norm or Angular Velocity Norm [deg/s]')
        title('Angular Velocity vs. Orbit, Bang-Bang Law')
        grid minor
        semilogy(tspan./T, 3*ones(1,length(tspan)),'k');
    end
    
end
    