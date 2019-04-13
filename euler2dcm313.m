function [dcm] = euler2dcm313(euler)
%EULER2DCM313 Convert 3-1-3 euler angles to a DCM
%   Detailed explanation goes here

psi = euler(1);
theta = euler(2);
phi = euler(3);

dcm = [cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(psi), cos(phi)*sin(psi) + sin(phi)*cos(theta)*cos(psi), sin(phi)*sin(theta);
       -sin(phi)*cos(psi) - cos(phi)*cos(theta)*sin(psi), -sin(phi)*sin(psi) + cos(phi)*cos(theta)*cos(psi), cos(phi)*sin(theta);
       sin(theta)*sin(psi), -sin(theta)*cos(psi), cos(theta)];
   
end

