function out_m = loveraTrackControl1(Bb, sigma, omega, control, Gamma)
%TRACKCONTROL Attitude Tracking control

global m_max;
global I;

sigma_RN = control.sigma_RN;
K_sigma = control.K_sigma;
%K_omega = control.K_omega;
epsilon = control.epsilon;
sigma_BR = sigma - sigma_RN;

K_omega = sqrt(K_sigma*min(eig(Gamma))^2/min(diag(I))*sqrt(cond(Gamma))) + 0.1;

S = tilde(Bb);

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