function out_m = loveraTrackControl1(Bb, sigma, omega, control)
%TRACKCONTROL Attitude Tracking control

global m_max;
global I;

sigma_RN = control.sigma_RN;
K_sigma = control.K_sigma;
K_omega = control.K_omega;
epsilon = control.epsilon;
sigma_BR = sigma - sigma_RN;

bhat = Bb./norm(Bb);

Gamma = (eye(3) - bhat*bhat');

S = [0, Bb(3), -Bb(2); -Bb(3), 0, Bb(1); Bb(2), -Bb(1), 0];

q = 2*sigma_BR/(1-norm(sigma_BR)^2);

u = -(epsilon^2*K_sigma*q + epsilon*K_omega*I*omega);

out_m = 1/norm(Bb)^2*S'*u;

% apply dipole constraint
for i = 1:3
    if abs(out_m(i)) > m_max
        out_m(i) = sign(out_m(i))*m_max;
    end
end

end