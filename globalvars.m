%GLOBALVARS
% Script which contains the global variables of the project

% List of globals
global tstep;
global I;
global beta_m_0; 
global M;
global gamma_m;
global omega_e;
global R_Earth;
global h;
global r;
global mu;
global Omega;
global arg_per;
global theta0;
global thetadot;
global T;
global inc;
global ecc;
global gain;
global m_max;
global use_gravity_gradient;
global R_c_o;

% Initialize globals
tstep = 0.1; % integration time-step
I = diag([1/12*4.8*(0.1^2+0.1^2), 1/12*4.8*(0.3^2+0.1^2), 1/12*4.8*(0.3^2+0.1^2)]); % [kg-m^2] inertia tensor
%I = diag([0.1,0.2,0.05]);
beta_m_0 = 0; % [rad]
M = 7.838e6; % [T-km^3] dipole moment
gamma_m = deg2rad(17); % [rad] tilt angle
omega_e = 7.2921159e-5; % [rad/s] earth rotation rate about n3
R_Earth= 6378; % [km] Mean earth radius
h = 405; % [km] altitude
r = (h + R_Earth); % [km] orbital radius
mu = 398600; % [km^2/s^2] gravitational parameter
Omega = 0; % [rad] right ascension of the ascending node
arg_per = 0;
inc = deg2rad(51.6388); % [rad] inclination
ecc = 0.0003449;
theta0 = deg2rad(80); % [rad] true anomaly (as argument of periapsis is not defined)
thetadot = sqrt(mu/r^3); % [rad/s] change in true anomaly
T = 2*pi*sqrt(r^3/mu); % [s] orbital period
gain = -1; % Default gain given by function
m_max = 0.11; % [A-m^2] maximum magnetic dipole from each rod (air torquer)
use_gravity_gradient = 1;
R_c_o = [0, 0, r]'; % orbital frame spacecraft position vector