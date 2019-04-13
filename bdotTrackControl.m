function out_m = bdotTrackControl(Bb, sigma, omega, control, beta_m)
%TRACKCONTROL Attitude Tracking control

% Global variables
global I;
global thetadot;
global inc;
global gamma_m;
global Omega;
global gain;
global m_max;

sigma_RN = control.sigma_RN;
K_sigma = control.K_sigma;
sigma_BR = sigma - sigma_RN;

bhat = Bb./norm(Bb);

    
xi_m = acos(cos(inc)*cos(gamma_m) + sin(inc)*sin(gamma_m)*cos(Omega - beta_m));

% Calculate gain if input == -1
if (gain == -1)
    K_omega = 2*thetadot*(1+sin(xi_m))*min(diag(I));
else
    K_omega = gain;
end

out_m = cross((eye(3) - bhat*bhat')*(K_omega*omega + K_sigma*sigma_BR), bhat/norm(Bb));

% apply dipole constraint
for i = 1:3
    if abs(out_m(i)) > m_max
        out_m(i) = sign(out_m(i))*m_max;
    end
end

end