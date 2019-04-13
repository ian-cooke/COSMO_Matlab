function out_m = trackControl(Bb, sigma, omega, control)
%TRACKCONTROL Attitude Tracking control

global m_max;

sigma_RN = control.sigma_RN;
K_sigma = control.K_sigma;
K_omega = control.K_omega;
sigma_BR = sigma - sigma_RN;

out_m = 1/norm(Bb)^2*cross(-K_sigma*sigma_BR + K_omega*omega, Bb);

% apply dipole constraint
for i = 1:3
    if abs(out_m(i)) > m_max
        out_m(i) = sign(out_m(i))*m_max;
    end
end

end

