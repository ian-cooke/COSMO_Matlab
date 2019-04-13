%MS4
% Milestone 4 Script
% Control Law Implementation

% Global Vars
global T;
global tstep;
global gain;

% Implement the two B-dot control laws and integrate the model
% Detumbling needs to occur within 3 orbits.
% Plot the magnetic dipole m for one of the initial conditions
% Case 1 - Init Condit 1 and Modulating b-dot control law
sigma0 = [0.3, 0.2, 0.4]'; % [none]
omega0 = deg2rad([15, 8, 12]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'mod';
[Xcase1, uCase1, mCase1] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase1 = Xcase1(1:3, :);
omegaCase1 = Xcase1(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase1 = findThresh(Xcase1(4:6, :));

f6 = figure(6);
plotState(f6, tspan, sigmaCase1, omegaCase1);

f7 = figure(7);
plotControl(f7, tspan, uCase1, mCase1);


% Case 2 - Init Condit 1 and bang-bang b-dot control law
sigma0 = [0.3, 0.2, 0.4]'; % [none]
omega0 = deg2rad([15, 8, 12]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'bang';
[Xcase2, uCase2, mCase2] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase2 = Xcase2(1:3, :);
omegaCase2 = Xcase2(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase2 = findThresh(Xcase2(4:6, :));

f8 = figure(8);
plotState(f8, tspan, sigmaCase2, omegaCase2);
f9 = figure(9);
plotControl(f9, tspan, uCase2, mCase2);

% Create difference in angular velocity components
f50 = figure(50);
set(f50, 'defaultaxesfontsize', 16)
xlabel('Time t [s]')
ylabel('Angular Velocity \omega [deg/s]')
hold on
grid minor
delomega = omegaCase2 - omegaCase1;
plot(tspan, delomega(1, :), 'Linewidth', 1)
plot(tspan, delomega(2, :), 'Linewidth', 1)
plot(tspan, delomega(3, :), 'Linewidth', 1)
leg = legend('\delta\omega_1', '\delta\omega_2', '\delta\omega_3');
set(leg,'location','best')

% Case 3 - Init Condit 2 and Modulating b-dot control law
sigma0 = [0.1, 0.1, 0.4]'; % [none]
omega0 = deg2rad([1, 12, 1]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'mod';
[Xcase3, uCase3, mCase3] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase3 = Xcase3(1:3, :);
omegaCase3 = Xcase3(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase3 = findThresh(Xcase3(4:6, :));

f10 = figure(10);
plotState(f10, tspan, sigmaCase3, omegaCase3);

f11 = figure(11);
plotControl(f11, tspan, uCase3, mCase3);

% Employ different gains
gain = 0.0005;
[XcaseGain005,uCaseGain005,mCaseGain005] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCaseGain005 = XcaseGain005(1:3, :);
omegaCaseGain005 = XcaseGain005(4:6, :);

f40 = figure(40);
plotState(f40, tspan, sigmaCaseGain005, omegaCaseGain005);

f41 = figure(41);
plotControl(f41, tspan, uCaseGain005, mCaseGain005);

gain = 0.001;
[XcaseGain01,uCaseGain01,mCaseGain01] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCaseGain01 = XcaseGain01(1:3, :);
omegaCaseGain01 = XcaseGain01(4:6, :);

f42 = figure(42);
plotState(f42, tspan, sigmaCaseGain01, omegaCaseGain01);

f43 = figure(43);
plotControl(f43, tspan, uCaseGain01, mCaseGain01);

% Case 4 - Init Condit 1 and bang-bang b-dot control law
gain = -1; % reset to default
sigma0 = [0.1, 0.1, 0.4]'; % [none]
omega0 = deg2rad([1, 12, 1]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'bang';
[Xcase4, uCase4, mCase4] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase4 = Xcase4(1:3, :);
omegaCase4 = Xcase4(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase4 = findThresh(Xcase4(4:6, :));

f12 = figure(12);
plotState(f12, tspan, sigmaCase4, omegaCase4);
f13 = figure(13);
plotControl(f13, tspan, uCase4, mCase4);

% Case 5 - Init Condit 3 and Modulating b-dot control law
sigma0 = [0.3, 0.2, 0.15]'; % [none]
omega0 = deg2rad([6, 4, 13]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'mod';
[Xcase5, uCase5, mCase5] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase5 = Xcase5(1:3, :);
omegaCase5 = Xcase5(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase5 = findThresh(Xcase5(4:6, :));

f14 = figure(14);
plotState(f14, tspan, sigmaCase5, omegaCase5);

f15 = figure(15);
plotControl(f15, tspan, uCase5, mCase5);


% Case 6 - Init Condit 1 and bang-bang b-dot control law
sigma0 = [0.3, 0.2, 0.15]'; % [none]
omega0 = deg2rad([6, 4, 13]'); % [rad/s]
x0 = [sigma0;omega0];
tspan = 0:tstep:3*T;
control = 'bang';
[Xcase6, uCase6, mCase6] = intmrprk4(tspan, I, x0, control, Hb);
sigmaCase6 = Xcase6(1:3, :);
omegaCase6 = Xcase6(4:6, :);

% find where angular velocity meets 3 deg/s threshold
threshCase6 = findThresh(Xcase6(4:6, :));

f16 = figure(16);
plotState(f16, tspan, sigmaCase6, omegaCase6);
f17 = figure(17);
plotControl(f17, tspan, uCase6, mCase6);


