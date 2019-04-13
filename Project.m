% 
% %% Milestone 1: Coordinate frame relations
% % Coordinate Frames:
% %N = (O, n1, n2, n3);
% %B = (H, b1, b2, b3);
% %H = (O, ir, itheta, ih);
% %M = (O, m1, m2, m3);
% %See Project Assignment for definition of coordinate frames.
% %Derivation of [BN] using (3-1-3) Euler Angles
% %B frame - (3-1-3) = (, , )
% % syms psi theta phi
% %     BN = euler2dcm313([psi, theta, phi])
% %Derivation of [HN] using (3-1-3) Euler Angles = (, , ).
% % syms Omega i nu
% % HN = euler2dcm313([0, i, nu])
% %Derivation of [MN] using (3-1-3) Euler Angles
% %M-frame - (3-1-3) = (, , 0) = (, , 0)
% % syms beta_m gamma_mm
% % MN = euler2dcm313([beta_m, gamma_mm, 0])
% %Derivation of [HM] using (3-1-3) Euler Angles
% % HM = HN*MN'
% % % Clear the variables
% % BN = [];
% % HN = [];
% 
% %% Milestone 2: Numerical Integrator
% % RK4 Numerical Integrator
% % State: X = [sigma, omega]
% % This function can be found at the very bottom of this file.
% 
% % Consistent tstep function
% global tstep;
% tstep = 0.1;
% 
% % Integrate the state with my integrator for 100 seconds
% sigma0 = [0.3, 0.2, 0.4]';
% omega0 = deg2rad([15, 8, 12]'); % [rad/s]
% x0 = [sigma0;omega0];
% u = 'none';
% global I;
% I = diag([3.5, 5, 8]); % [kg-m^2]
% tspan = 0:tstep:100; % [s]
% 
% [state, ~, ~] = intmrprk4(tspan, I, x0, u, 0);
% sigma = state(1:3, :);
% omega = state(4:6, :);
% 
% f1 = figure(1);
% set(f1, 'defaultaxesfontsize', 16)
% set(f1, 'Visible', 'on')
% subplot(2,1,1)
% xlabel('Time t [s]')
% ylabel('MRP \sigma')
% hold on
% grid minor
% plot(tspan, sigma(1, :), 'Linewidth', 3)
% plot(tspan, sigma(2, :), 'Linewidth', 3)
% plot(tspan, sigma(3, :), 'Linewidth', 3)
% leg = legend('\sigma_1', '\sigma_2', '\sigma_3');
% set(leg, 'Location', 'best')
% hold off
% 
% subplot(2,1,2)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega')
% hold on
% grid minor
% plot(tspan, omega(1, :), 'Linewidth', 3)
% plot(tspan, omega(2, :), 'Linewidth', 3)
% plot(tspan, omega(3, :), 'Linewidth', 3)
% leg = legend('\omega_1', '\omega_2', '\omega_3');
% set(leg, 'Location', 'best')
% hold off
% % Plot Angular Momentum and Kinetic Energy (should be conserved)
% 
% H = zeros(1, length(tspan));
% KE = H;
% for i = 1:length(tspan)
%     Iomeg = I*omega(:, i);
%     H(i) = norm(Iomeg);
%     KE(i) = 1/2*norm(dot(Iomeg, omega(:, i)));
% end
% 
% f3 = figure(3);
% set(f3, 'defaultaxesfontsize', 16)
% set(f3, 'Visible', 'on')
% xlabel('Time t [s]')
% ylabel('Angular Momentum H [kg-m^2/s]')
% hold on
% yyaxis left
% plot(tspan, H, 'Linewidth', 3)
% ylim([2.0 2.1])
% yyaxis right
% plot(tspan, KE, 'Linewidth', 3)
% ylabel('Kinetic Energy T [J]')
% ylim([0.3 0.4])
% hold off
% % The angular momentum magnitude and kinetic energy are indeed constant.
% 
% %% Milestone 3: Geomagnetic Modeling
% % Modeling the tilted magnetic dipole using the coordinate frame representations
% global beta_m_0; 
% global M;
% global gamma_m;
% global omega_e;
% beta_m_0 = 0; % [rad]
% M = 7.838e6; % [T-km^3] dipole moment
% gamma_m = deg2rad(17); % [rad] tilt angle
% omega_e = 7.2921159e-5; % [rad/s] earth rotation rate about n3
% 
% % Orbital parameters
% global R_Earth;
% global h;
% global r;
% global mu;
% R_Earth= 6378; % [km] Mean earth radius
% h = 450; % [km] altitude
% r = (h + R_Earth); % [m] orbital radius
% mu = 398600; % [km^2/s^2] gravitational parameter
% 
% % LVLH reference frame angles
% global Omega;
% global theta0;
% global thetadot;
% global T;
% global inc;
% Omega = 0; % [rad] right ascension of the ascending node
% inc = deg2rad(45); % [rad] inclination
% theta0 = 0; % [rad] true anomaly (as argument of periapsis is not defined)
% thetadot = sqrt(mu/r^3); % [rad/s] change in true anomaly
% T = 2*pi*sqrt(r^3/mu); % [s] orbital period
% 
% % Plot the magnetic field  as seen by the orbit frame H over one orbit
% % No rk4 needed - no differential equation - simply nee a loop
% tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
% % Propagate the state
% Hb = calcMagField(tspan);
% 
% f4 = figure(4);
% set(f4, 'Visible', 'on')
% set(f4, 'defaultaxesfontsize', 16)
% hold on
% grid minor
% xlabel('Time t [s]')
% xlim([0 T])
% ylabel('Magnetic Field Strength [T]')
% plot(tspan, Hb(1, :), 'Linewidth', 3)
% plot(tspan, Hb(2, :), 'Linewidth', 3)
% plot(tspan, Hb(3, :), 'Linewidth', 3)
% leg = legend('^Hb_1', '^Hb_2', '^Hb_3');
% set(leg, 'Location', 'best')
% 
% %Plot the Magnetic Field as seen by the body frame for 100 seconds
% % Just need to convert H coordinates to B coordinates
% % Bb = [BN][HN]'*Hb
% tspan = 0:tstep:100;
% h = tspan(2) - tspan(1);
% t = tspan(1);
% ind = 1;
% while (t <= tspan(end))
%     % [BN] DCM using MRPs calculated in Milestone 2
%     s = sigma(:, ind);
%     sq = s'*s;
%     BN = eye(3) + (8*tilde(s)*tilde(s) - 4*(1-sq)*tilde(s))/(1+sq)^2;
%                  
%     % [HN] DCM
%     theta = theta0 + thetadot*t;
%     HN = [cos(Omega)*cos(theta) - sin(Omega)*cos(inc)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(inc)*sin(theta), sin(inc)*sin(theta);
%           -cos(Omega)*sin(theta) - sin(Omega)*cos(inc)*cos(theta), cos(Omega)*cos(inc)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(inc);
%           sin(Omega)*sin(inc), -cos(Omega)*sin(inc), cos(inc)];
%       
%     % [BH] DCM
%     BH = BN*HN';
%     Bb(:, ind) = BH*Hb(:, ind);
%     
% t = t + h;
% ind = ind + 1;
% end
% 
% % Plot it
% f5 = figure(5);
% set(f5, 'Visible', 'on')
% set(f5, 'defaultaxesfontsize', 16)
% hold on
% grid minor
% xlabel('Time t [s]')
% ylabel('Magnetic Field Strength [T]')
% plot(tspan, Bb(1, :), 'Linewidth', 3)
% plot(tspan, Bb(2, :), 'Linewidth', 3)
% plot(tspan, Bb(3, :), 'Linewidth', 3)
% leg = legend('^Bb_1', '^Bb_2', '^Bb_3');
% set(leg, 'Location', 'best')
% 
% %% Milestone 4: Control Law Implementation
% % Implement the two B-dot control laws and integrate the model
% % Detumbling needs to occur within 3 orbits.
% % Use the following initial condtions for each control law
% % deg/s
% % deg/s
% % deg/s
% % Plot the magnetic dipole m for one of the initial conditions
% % Case 1 - Init Condit 1 and Modulating b-dot control law
% sigma0 = [0.3, 0.2, 0.4]'; % [none]
% omega0 = deg2rad([15, 8, 12]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'mod';
% [Xcase1, uCase1, mCase1] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase1 = Xcase1(1:3, :);
% omegaCase1 = Xcase1(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase1 = findThresh(Xcase1(4:6, :));
% 
% f6 = figure(6);
% plotState(f6, tspan, sigmaCase1, omegaCase1);
% 
% f7 = figure(7);
% plotControl(f7, tspan, uCase1, mCase1);
% 
% 
% % Case 2 - Init Condit 1 and bang-bang b-dot control law
% sigma0 = [0.3, 0.2, 0.4]'; % [none]
% omega0 = deg2rad([15, 8, 12]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'bang';
% [Xcase2, uCase2, mCase2] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase2 = Xcase2(1:3, :);
% omegaCase2 = Xcase2(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase2 = findThresh(Xcase2(4:6, :));
% 
% f8 = figure(8);
% plotState(f8, tspan, sigmaCase2, omegaCase2);
% f9 = figure(9);
% plotControl(f9, tspan, uCase2, mCase2);
% 
% % Case 3 - Init Condit 2 and Modulating b-dot control law
% sigma0 = [0.1, 0.1, 0.4]'; % [none]
% omega0 = deg2rad([1, 12, 1]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'mod';
% [Xcase3, uCase3, mCase3] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase3 = Xcase3(1:3, :);
% omegaCase3 = Xcase3(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase3 = findThresh(Xcase3(4:6, :));
% 
% f10 = figure(10);
% plotState(f10, tspan, sigmaCase3, omegaCase3);
% 
% f11 = figure(11);
% plotControl(f11, tspan, uCase3, mCase3);
% 
% 
% % Case 4 - Init Condit 1 and bang-bang b-dot control law
% sigma0 = [0.1, 0.1, 0.4]'; % [none]
% omega0 = deg2rad([1, 12, 1]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'bang';
% [Xcase4, uCase4, mCase4] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase4 = Xcase4(1:3, :);
% omegaCase4 = Xcase4(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase4 = findThresh(Xcase4(4:6, :));
% 
% f12 = figure(12);
% plotState(f12, tspan, sigmaCase4, omegaCase4);
% f13 = figure(13);
% plotControl(f13, tspan, uCase4, mCase4);
% 
% % Case 5 - Init Condit 3 and Modulating b-dot control law
% sigma0 = [0.3, 0.2, 0.15]'; % [none]
% omega0 = deg2rad([6, 4, 13]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'mod';
% [Xcase5, uCase5, mCase5] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase5 = Xcase5(1:3, :);
% omegaCase5 = Xcase5(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase5 = findThresh(Xcase5(4:6, :));
% 
% f14 = figure(14);
% plotState(f14, tspan, sigmaCase5, omegaCase5);
% 
% f15 = figure(15);
% plotControl(f15, tspan, uCase5, mCase5);
% 
% 
% % Case 6 - Init Condit 1 and bang-bang b-dot control law
% sigma0 = [0.3, 0.2, 0.15]'; % [none]
% omega0 = deg2rad([6, 4, 13]'); % [rad/s]
% x0 = [sigma0;omega0];
% tspan = 0:tstep:3*T;
% control = 'bang';
% [Xcase6, uCase6, mCase6] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase6 = Xcase6(1:3, :);
% omegaCase6 = Xcase6(4:6, :);
% 
% % find where angular velocity meets 3 deg/s threshold
% threshCase6 = findThresh(Xcase6(4:6, :));
% 
% f16 = figure(16);
% plotState(f16, tspan, sigmaCase6, omegaCase6);
% f17 = figure(17);
% plotControl(f17, tspan, uCase6, mCase6);
% 
% %% Milestone 5: Orbit Inclination v/s control performance
% % Case 1 - 15 deg inclination and modulating control
% inc = deg2rad(15);
% % Re-calculate the magnetic field for new inclination
% tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
% Hb = calcMagField(tspan);
% 
% % Now integrate
% sigma0 = [0.3,0.2,0.4]'; % [none]
% omega0 = deg2rad([15, 8, 12]'); % [rad/s]
% x0 = [sigma0;omega0];
% control = 'mod';
% [Xcase7, uCase7, mCase7] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase7 = Xcase7(1:3, :);
% omegaCase7 = Xcase7(4:6, :);
% 
% % Only need to plot control error
% f18 = figure(18);
% set(f18, 'defaultaxesfontsize', 16)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega [rad/s]')
% hold on
% grid minor
% plot(tspan, omegaCase7(1, :), 'Linewidth', 2)
% plot(tspan, omegaCase7(2, :), 'Linewidth', 2)
% plot(tspan, omegaCase7(3, :), 'Linewidth', 2)
% % threshold
% plot(tspan, deg2rad(3)*ones(1,length(tspan)),'k')
% plot(tspan, -deg2rad(3)*ones(1,length(tspan)),'k')
% leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
% set(leg, 'Location', 'best')
% 
% % Case 2 - 15 deg inclination and bang-bang control
% control = 'bang';
% [Xcase8, uCase8, mCase8] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase8 = Xcase8(1:3, :);
% omegaCase8 = Xcase8(4:6, :);
% 
% % Only need to plot control error
% f19 = figure(19);
% set(f19, 'defaultaxesfontsize', 16)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega [rad/s]')
% hold on
% grid minor
% plot(tspan, omegaCase8(1, :), 'Linewidth', 2)
% plot(tspan, omegaCase8(2, :), 'Linewidth', 2)
% plot(tspan, omegaCase8(3, :), 'Linewidth', 2)
% % threshold
% plot(tspan, deg2rad(3)*ones(1,length(tspan)),'k')
% plot(tspan, -deg2rad(3)*ones(1,length(tspan)),'k')
% leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
% set(leg, 'Location', 'best')
% 
% % Case 3 - 105 deg inclination with modulating control
% inc = deg2rad(105);
% % Re-calculate the magnetic field for new inclination
% tspan = 0:tstep:3*T; % [s] time propagation vector - 3 orbits to allow adequate time for detumbling
% Hb = calcMagField(tspan);
% 
% % Now integrate
% control = 'mod';
% [Xcase9, uCase9, mCase9] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase9 = Xcase9(1:3, :);
% omegaCase9 = Xcase9(4:6, :);
% 
% % Only need to plot control error
% f20 = figure(20);
% set(f20, 'defaultaxesfontsize', 16)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega [rad/s]')
% hold on
% grid minor
% plot(tspan, omegaCase9(1, :), 'Linewidth', 2)
% plot(tspan, omegaCase9(2, :), 'Linewidth', 2)
% plot(tspan, omegaCase9(3, :), 'Linewidth', 2)
% % threshold
% plot(tspan, deg2rad(3)*ones(1,length(tspan)),'k')
% plot(tspan, -deg2rad(3)*ones(1,length(tspan)),'k')
% leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
% set(leg, 'Location', 'best')
% 
% % Case 4 - 105 deg inclination with bang-bang control
% control = 'bang';
% [Xcase10, uCase10, mCase10] = intmrprk4(tspan, I, x0, control, Hb);
% sigmaCase10 = Xcase10(1:3, :);
% omegaCase10 = Xcase10(4:6, :);
% 
% % Only need to plot control error
% f21 = figure(21);
% set(f21, 'defaultaxesfontsize', 16)
% xlabel('Time t [s]')
% ylabel('Angular Velocity \omega [rad/s]')
% hold on
% grid minor
% plot(tspan, omegaCase10(1, :), 'Linewidth', 2)
% plot(tspan, omegaCase10(2, :), 'Linewidth', 2)
% plot(tspan, omegaCase10(3, :), 'Linewidth', 2)
% % threshold
% plot(tspan, deg2rad(3)*ones(1,length(tspan)),'k')
% plot(tspan, -deg2rad(3)*ones(1,length(tspan)),'k')
% leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
% set(leg, 'Location', 'best')
% %% Milestone 6
% 
% 
% %% Script Functions
% function [X, u, m] = intmrprk4(tspan, I, x0, control, Hb)
% %INTMRPRK4 Numerically integrate the Modified Rodriguez Parameters and
% %angular velocity vector using the rigid body equations of motion and the 
% %RK4 numerical method. 
% % 
% % Purpose
% %   To numerically integrate the state x = [sigma, omega]'
% %
% % Inputs
% %   tspan - vector of time points, evenly spaced
% %   I - inertia tensor of rigid body (aligned with princip. axis frame)
% %   x0 - initial state, [sigma0, omega0] (must be a column vector)
% %   u - control law
% % Outputs
% %   X - matrix of states propagated with t
% %
% % Author(s):
% %   Ian Coooke
% %
% % Created
% %   8 Apr 2018
% % Modified
% %   23 Apr 2018
% % Log
% %   8 Apr 2018
% %   23 Apr 2018
% %       Added control
% %       
% %-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-%
% 
% % Define global variables
% global theta0;
% global thetadot;
% global Omega;
% global inc;
% global beta_m_0;
% global omega_e;
% 
% % Extract the states
% sigma0 = x0(1:3);
% omega0 = x0(4:6);
% 
% % Integration setup
% n = length(tspan);
% 
% sigma = zeros(3, n);
% omega = sigma;
% 
% u = sigma;
% m = sigma;
% 
% Bb = sigma;
% Bb(:, end+1) = [0,0,0]';
% 
% sigma(:, 1) = sigma0;
% omega(:, 1) = omega0;
% 
% h = tspan(2) - tspan(1);
% t = tspan(1);
% 
% % Integration loop
% tic
% for ind=1:n-1
%     
%     % current states
%     s = sigma(:, ind);
%     w = omega(:, ind);
%     
%     % MRPs
%     % k1
%     sq = s'*s; %sigma^2
%     Q = tilde(s); % tilde(sigma)
%     B = (1 - sq)*eye(3) + 2*Q + 2*s*s'; % B-matrix
%     k1 = 1/4*B*w;
%     
%     % k2
%     s2 = s + h*k1/2;
%     sq = s2'*s2;
%     Q = tilde(s2);
%     B = (1- sq)*eye(3) + 2*Q + 2*s2*s2';
%     k2 = 1/4*B*w;
%     
%     % k3
%     s3 = s + h*k2/2;
%     sq = s3'*s3;
%     Q = tilde(s3);
%     B = (1- sq)*eye(3) + 2*Q + 2*s3*s3';
%     k3 = 1/4*B*w;
%     
%     % k4
%     s4 = s + h*k3;
%     sq = s4'*s4;
%     Q = tilde(s4);
%     B = (1- sq)*eye(3) + 2*Q + 2*s4*s4';
%     k4 = 1/4*B*w;
%     
%     % next state
%     sigma(:, ind+1) = sigma(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
%     
%     % MRP switching
%     if (norm(sigma(:, ind+1)) > 1)
%         sigma(:, ind+1) = -sigma(:, ind+1)/norm(sigma(:, ind+1))^2;
%     end
%     
%     % control law
%     % only calculate magnetic field if you're using control
%     if strcmp(control, 'mod') || strcmp(control, 'bang')
%          s = sigma(:, ind);
%         % Magnetic field in B components (H components already calculated)
%         BN = eye(3) + (tilde(s)*(8*tilde(s) - 4*(1-s'*s)))/(1+s'*s)^2;
%                      
%         % [HN] DCM
%         theta = theta0 + thetadot*t;
%         beta_m = beta_m_0 + omega_e*t;
%         HN = [cos(Omega)*cos(theta) - sin(Omega)*cos(inc)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(inc)*sin(theta), sin(inc)*sin(theta);
%               -cos(Omega)*sin(theta) - sin(Omega)*cos(inc)*cos(theta), cos(Omega)*cos(inc)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(inc);
%               sin(Omega)*sin(inc), -cos(Omega)*sin(inc), cos(inc)];
%           
%         % [BH] DCM and derivative Bb'
%         BH = BN*HN';
%         Bb(:, ind) = BH*Hb(:, ind);
%         Bbprime = BH*diff([Hb(:, ind), Hb(:, ind+1)],1,2);
%         
%         % modulating b-dot control law
%         if strcmp(control, 'mod')
%             m(:, ind) = modControl(Bb(:, ind), beta_m, omega(:, ind));
%             
%         % bang-bang b-dot control law
%         elseif strcmp(control, 'bang')
%             % Numerical derivative
%             m(:, ind) = bangControl(Bbprime);
%         end
%     elseif strcmp(control, 'none')
%         m(:, ind) = [0,0,0]';
%         Bb(:, ind) = m(:, ind);
%     else
%         error('invalid control law')
%     end
%     
%     % Calculate control
%     u(:, ind) = cross(m(:, ind), Bb(:, ind)); 
%     
%     % omega
%     % k1
%     k1 = I^-1*(-tilde(w)*I*w + u(:, ind));
%     
%     % k2
%     w2 = w + h/2*k1;
%     k2 = I^-1*(-tilde(w2)*I*w2 + u(:, ind));
%     
%     % k3
%     w3 = w + h/2*k2;
%     k3 = I^-1*(-tilde(w3)*I*w3 + u(:, ind));
%     
%     % k4
%     w4 = w + h*k3;
%     k4 = I^-1*(-tilde(w4)*I*w4 + u(:, ind));
%     
%     % next state
%     omega(:, ind+1) = omega(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
%     
% t = t + h;
% end
% toc
% 
% X = [sigma; omega];
% 
% end
% 
% function [m] = modControl(Bb, beta_m, omega)
% %MODCONTROL Numerically calculate the magnetic dipole to command the control torque
% % 
% % Purpose
% %   To numerically calculate the magnetic dipole moment
% %
% % Inputs
% %   Bb - Magnetic field at instant in time as seen by satellite frame
% %   beta_m - angle between magnetic frame and inertial frame
% %   omega - current angular velocity vector as seen by inertial frame
% % Outputs
% %   m - magnetic dipole
% %
% % Author(s):
% %   Ian Coooke
% %
% % Created
% %   23 Apr 2018
% % Modified
% %   23 Apr 2018
% % Log
% %   23 Apr 2018
% %       
% %-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-%
%     % Global variables
%     global I;
%     global thetadot;
%     global inc;
%     global gamma_m;
%     global Omega;
%     
%     m_max = 3;
%     
%     bhat = Bb./norm(Bb);
%     
%     xi_m = acos(cos(inc)*cos(gamma_m) + sin(inc)*sin(gamma_m)*cos(Omega - beta_m));
%     k_w = 2*thetadot*(1+sin(xi_m))*min(diag(I));
%     %k_w = 0.005;
%     m = cross(-k_w/norm(Bb)*bhat, (eye(3) - Bb*Bb')*omega);
%     
%     % Enforce dipole saturation
%     for j = 1:3
%         if abs(m(j)) > m_max
%             m(j) = sign(m(j))*m_max;
%         end
%     end
% end
% 
% function [m] = bangControl(bp)
%      m_max = 3;
%      m = m_max*sign(bp);
% end
% 
% function [res] = findThresh(omega)
%     thresh = deg2rad(3);
%     res = zeros(3,1);
%     for j = 1:length(omega)
%         if omega(1,j) < thresh
%             res(1) = j;
%         elseif omega(2,j) < thresh
%             res(2) = j;
%         elseif omega(3,j) < thresh
%             res(3) = j;
%         end
%         res = max(res);
%     end
% end
% 
% function plotState(fig, tspan, sigma, omega)
%     set(fig, 'Visible', 'on')
%     set(fig, 'defaultaxesfontsize', 16)
% 
%     subplot(2,1,1)
%     xlabel('Time t [s]')
%     ylabel('MRP \sigma')
%     hold on
%     grid minor
%     plot(tspan, sigma(1, :), 'Linewidth', 2)
%     plot(tspan, sigma(2, :), 'Linewidth', 2)
%     plot(tspan, sigma(3, :), 'Linewidth', 2)
%     leg = legend('\sigma_1', '\sigma_2', '\sigma_3');
%     set(leg, 'Location', 'best')
% 
%     subplot(2,1,2)
%     xlabel('Time t [s]')
%     ylabel('Angular Velocity \omega [rad/s]')
%     hold on
%     grid minor
%     plot(tspan, omega(1, :), 'Linewidth', 2)
%     plot(tspan, omega(2, :), 'Linewidth', 2)
%     plot(tspan, omega(3, :), 'Linewidth', 2)
%     % threshold
%     plot(tspan, deg2rad(3)*ones(1,length(tspan)),'k')
%     plot(tspan, -deg2rad(3)*ones(1,length(tspan)),'k')
%     leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold');
%     set(leg, 'Location', 'best')
% end
% 
% function plotControl(fig, tspan, u, m)
%     set(fig, 'Visible', 'on')
%     set(fig, 'defaultaxesfontsize', 16)
% 
%     subplot(2,1,1)
%     xlabel('Time t [s]')
%     ylabel('Control Torque u [N-m]')
%     hold on
%     grid minor
%     plot(tspan, u(1, :), 'Linewidth', 2)
%     plot(tspan, u(2, :), 'Linewidth', 2)
%     plot(tspan, u(3, :), 'Linewidth', 2)
%     leg = legend('u_1', 'u_2', 'u_3');
%     set(leg, 'Location', 'best')
% 
%     subplot(2,1,2)
%     xlabel('Time t [s]')
%     ylabel('Magnetic Dipole m [T]')
%     hold on
%     grid minor
%     plot(tspan, m(1, :), 'Linewidth', 2)
%     plot(tspan, m(2, :), 'Linewidth', 2)
%     plot(tspan, m(3, :), 'Linewidth', 2)
%     leg = legend('m_1', 'm_2', 'm_3');
%     set(leg, 'Location', 'best')
% end
% 
% function [m] = tilde(v)
% %TILDE Skew symmetric matrix operator
% 
% m = [0, -v(3), v(2);
%      v(3), 0, -v(1);
%      -v(2), v(1), 0];
% 
% end
% 
% function [Hb] = calcMagField(tspan)
%     % Global Variables
%     global inc;
%     global beta_m_0;
%     global omega_e;
%     global thetadot;
%     global gamma_m;
%     global Omega;
%     global h;
%     global M;
%     global r;
%     
%     % Initialize vector
%     Hb = zeros(3, length(tspan));
%     % initial time
%     t = tspan(1);
%     % Propagate the state
%     ind = 1;
%     while (ind < length(tspan)+1)
%         % Changing params
%         beta_m = beta_m_0 + omega_e*t;
%         xci_m = acos(cos(inc)*cos(gamma_m) + sin(inc)*sin(gamma_m)*cos(Omega - beta_m));
%         eta_m = asin(sin(gamma_m)*sin(Omega - beta_m)/sin(xci_m));
% 
%         % Magnetic field
%         Hb(:, ind) = M/r^3*[cos(thetadot*t - eta_m)*sin(xci_m);
%                             cos(xci_m);
%                             -2*sin(thetadot*t - eta_m)*sin(xci_m)];
% 
%     ind = ind + 1;
%     t = t + h;
%     end
%     
% end