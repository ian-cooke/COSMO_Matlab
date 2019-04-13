%MS1
% Milestone 1 Script
% Coordinate Frame Relations

% Coordinate Frames:
%N = (O, n1, n2, n3);
%B = (H, b1, b2, b3);
%H = (O, ir, itheta, ih);
%M = (O, m1, m2, m3);
%See Project Assignment for definition of coordinate frames.
%Derivation of [BN] using (3-1-3) Euler Angles
%B frame - (3-1-3) = (, , )
% syms psi theta phi
%     BN = euler2dcm313([psi, theta, phi])
%Derivation of [HN] using (3-1-3) Euler Angles = (, , ).
% syms Omega i nu
% HN = euler2dcm313([0, i, nu])
%Derivation of [MN] using (3-1-3) Euler Angles
%M-frame - (3-1-3) = (, , 0) = (, , 0)
% syms beta_m gamma_mm
% MN = euler2dcm313([beta_m, gamma_mm, 0])
%Derivation of [HM] using (3-1-3) Euler Angles
% HM = HN*MN'
% % Clear the variables
% BN = [];
% HN = [];