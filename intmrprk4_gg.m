function [X] = intmrprk4_gg(tspan, I, x0)
%INTMRPRK4_GG Numerically integrate the Modified Rodriguez Parameters and
%angular velocity vector using the rigid body equations of motion and the 
%RK4 numerical method in the presence of a gravity gradient.
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
global thetadot;
global tstep;

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
    k1 = 1/4*B*w - thetadot/4*[2*(s(1)*s(2) + s(3)); 2*s(2)^2 + 1 - sq; 2*(s(2)*s(3) - s(1))];
    
    % k2
    s2 = s + h*k1/2;
    sq = s2'*s2;
    Q = tilde(s2);
    B = (1- sq)*eye(3) + 2*Q + 2*s2*s2';
    k2 = 1/4*B*w - thetadot/4*[2*(s2(1)*s2(2) + s2(3)); 2*s2(2)^2 + 1 - sq; 2*(s2(2)*s2(3) - s2(1))];
    
    % k3
    s3 = s + h*k2/2;
    sq = s3'*s3;
    Q = tilde(s3);
    B = (1- sq)*eye(3) + 2*Q + 2*s3*s3';
    k3 = 1/4*B*w - thetadot/4*[2*(s3(1)*s3(2) + s3(3)); 2*s3(2)^2 + 1 - sq; 2*(s3(2)*s3(3) - s3(1))];
    
    % k4
    s4 = s + h*k3;
    sq = s4'*s4;
    Q = tilde(s4);
    B = (1- sq)*eye(3) + 2*Q + 2*s4*s4';
    k4 = 1/4*B*w - thetadot/4*[2*(s4(1)*s4(2) + s4(3)); 2*s4(2)^2 + 1 - sq; 2*(s4(2)*s4(3) - s4(1))];
    
    % next state
    sigma(:, ind+1) = sigma(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % MRP switching
    if (norm(sigma(:, ind+1)) > 1)
        sigma(:, ind+1) = -sigma(:, ind+1)/norm(sigma(:, ind+1))^2;
    end
    
    % omega
    % k1
    k1 = I^-1*(-tilde(w)*I*w);
    
    % k2
    w2 = w + h/2*k1;
    k2 = I^-1*(-tilde(w2)*I*w2);
    
    % k3
    w3 = w + h/2*k2;
    k3 = I^-1*(-tilde(w3)*I*w3);
    
    % k4
    w4 = w + h*k3;
    k4 = I^-1*(-tilde(w4)*I*w4);
    
    % next state
    omega(:, ind+1) = omega(:, ind) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
t = t + h;
end
toc

X = [sigma; omega];

end