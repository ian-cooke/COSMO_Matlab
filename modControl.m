function [m] = modControl(Bb, beta_m, omega)
%MODCONTROL Numerically calculate the magnetic dipole to command the control torque
% 
% Purpose
%   To numerically calculate the magnetic dipole moment
%
% Inputs
%   Bb - Magnetic field at instant in time as seen by satellite frame
%   beta_m - angle between magnetic frame and inertial frame
%   omega - current angular velocity vector as seen by inertial frame
% Outputs
%   m - magnetic dipole
%
% Author(s):
%   Ian Coooke
%
% Created
%   23 Apr 2018
% Modified
%   23 Apr 2018
% Log
%   23 Apr 2018
%       
%-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-ASEN5010-%
    % Global variables
    global I;
    global thetadot;
    global inc;
    global gamma_m;
    global Omega;
    global gain;
    global m_max;
    
    bhat = Bb./norm(Bb);
    
    xi_m = acos(cos(inc)*cos(gamma_m) + sin(inc)*sin(gamma_m)*cos(Omega - beta_m));
    
    % Calculate gain if input == -1
    if (gain == -1)
        k_w = 2*thetadot*(1+sin(xi_m))*min(diag(I));
    else
        k_w = gain;
    end
    m = cross(-k_w/norm(Bb)*bhat, (eye(3) - bhat*bhat')*omega);
    
    % Enforce dipole saturation
    for j = 1:3
        if abs(m(j)) > m_max
            m(j) = sign(m(j))*m_max;
        end
    end
end