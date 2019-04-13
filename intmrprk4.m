function [X, u, m] = intmrprk4(tspan, I, x0, control, Hb)
%INTMRPRK4 Numerically integrate the Modified Rodriguez Parameters and
%angular velocity vector using the rigid body equations of motion and the 
%RK4 numerical method. 
% 
% Purpose
%   To numerically integrate the state x = [sigma, omega]'
%
% Inputs
%   tspan - vector of time points, evenly spaced
%   I - inertia tensor of rigid body (aligned with princip. axis frame)
%   x0 - initial state, [sigma0, omega0] (must be a column vector)
%   u - control law
% Outputs
%   X - matrix of states propagated with t
%
% Author(s):
%   Ian Coooke
%
% Created
%   8 Apr 2018
% Modified
%   23 Apr 2018
% Log
%   8 Apr 2018
%   23 Apr 2018
%       Added control
%       
%-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-%

% Define global variables
global theta0;
global thetadot;
global Omega;
global inc;
global beta_m_0;
global omega_e;
global tstep;
global use_gravity_gradient;
global mu;
global r;

% Extract the states
sigma0 = x0(1:3);
omega0 = x0(4:6);

% Integration setup
n = length(tspan);

sigma = zeros(3, n);
omega = sigma;

u = sigma;
m = sigma;

sigma(:, 1) = sigma0;
omega(:, 1) = omega0;

h = tstep;
t = tspan(1);

% Integration loop
tic
for ind=1:n-1
    
    % current states
    s = sigma(:, ind);
    w = omega(:, ind);
    
    % MRPs
    % k1
    sq = s'*s; %sigma^2
    Q = tilde(s); % tilde(sigma)
    B = (1 - sq)*eye(3) + 2*Q + 2*s*s'; % B-matrix
    k1 = 1/4*B*w;
    
    % k2
    s2 = s + h*k1/2;
    sq = s2'*s2;
    Q = tilde(s2);
    B = (1- sq)*eye(3) + 2*Q + 2*s2*s2';
    k2 = 1/4*B*w;
    
    % k3
    s3 = s + h*k2/2;
    sq = s3'*s3;
    Q = tilde(s3);
    B = (1- sq)*eye(3) + 2*Q + 2*s3*s3';
    k3 = 1/4*B*w;
    
    % k4
    s4 = s + h*k3;
    sq = s4'*s4;
    Q = tilde(s4);
    B = (1- sq)*eye(3) + 2*Q + 2*s4*s4';
    k4 = 1/4*B*w;
    
    % next state
    sigma(:, ind+1) = sigma(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % MRP switching
    if (norm(sigma(:, ind+1)) > 1)
        sigma(:, ind+1) = -sigma(:, ind+1)/norm(sigma(:, ind+1))^2;
    end
    
    % gravity gradient as an external torque
    if use_gravity_gradient
    end
    
    % control law
    % only calculate magnetic field if you're using control
    if strcmp(control.name, 'mod') || strcmp(control.name, 'bang') || strcmp(control.name, 'track') || strcmp(control.name, 'bdotTrack') || strcmp(control.name, 'loveraTrack1')
         s = sigma(:, ind);
        sq = s'*s;
        BN = eye(3) + (8*tilde(s)*tilde(s) - 4*(1-sq)*tilde(s))/(1+sq)^2;

        % [HN] DCM
        theta = theta0 + thetadot*t;
        beta_m = beta_m_0 + omega_e*t;
        HN = [cos(Omega)*cos(theta) - sin(Omega)*cos(inc)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(inc)*sin(theta), sin(inc)*sin(theta);
              -cos(Omega)*sin(theta) - sin(Omega)*cos(inc)*cos(theta), cos(Omega)*cos(inc)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(inc);
              sin(Omega)*sin(inc), -cos(Omega)*sin(inc), cos(inc)];

        % [BH] DCM
        BH = BN*HN';
        Bb = BH*Hb(:, ind);
        
        % modulating b-dot control law
        if strcmp(control.name, 'mod')
            m(:, ind) = modControl(Bb, beta_m, omega(:, ind));
            
        % bang-bang b-dot control law
        elseif strcmp(control.name, 'bang')
            % Numerical derivative
            Homega_HN = [0, 0, thetadot]'; % H-frame components
            Bomega_BH = omega(:, ind) - BH*Homega_HN;
            Bbprime = (-tilde(Bomega_BH)*BH)*Hb(:, ind)...
                + BH*diff([Hb(:, ind), Hb(:, ind+1)],1,2)/h;
            m(:, ind) = bangControl(Bbprime);
            
        % sigma track control
        elseif strcmp(control.name, 'track')
            m(:, ind) = trackControl(Bb, sigma(:, ind), omega(:, ind), control);
            
        % bdot track contorl    
        elseif strcmp(control.name, 'bdotTrack')
            m(:, ind) = bdotTrackControl(Bb, sigma(:, ind), omega(:, ind), control);
            
        % first control law from Lovera paper
        elseif strcmp(control.name, 'loveraTrack1')
            m(:, ind) = loveraTrackControl1(Bb, sigma(:, ind), omega(:, ind), control);
            
        % second control law from Lovera paper
        elseif strcmp(control.name, 'loveraTrack2')
            m(:, ind) = loveraTrackControl2(Bb, sigma(:, ind), omega(:, ind), control);
        end
        
        
    elseif strcmp(control.name, 'none')
        m(:, ind) = [0,0,0]';
        Bb = m(:, ind);
    else
        error('invalid control law')
    end
    
    % Calculate control
    u(:, ind) = cross(m(:, ind), Bb); 
    
    % omega
    % k1
    k1 = I^-1*(-tilde(w)*I*w + u(:, ind));
    
    % k2
    w2 = w + h/2*k1;
    k2 = I^-1*(-tilde(w2)*I*w2 + u(:, ind));
    
    % k3
    w3 = w + h/2*k2;
    k3 = I^-1*(-tilde(w3)*I*w3 + u(:, ind));
    
    % k4
    w4 = w + h*k3;
    k4 = I^-1*(-tilde(w4)*I*w4 + u(:, ind));
    
    % next state
    omega(:, ind+1) = omega(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
t = t + h;
end
toc

X = [sigma; omega];

end