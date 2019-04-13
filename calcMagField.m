function [Hb] = calcMagField(tspan)
    % Global Variables
    global inc;
    global beta_m_0;
    global omega_e;
    global thetadot;
    global gamma_m;
    global Omega;
    global M;
    global r;
    global tstep;
    
    % Initialize vector
    Hb = zeros(3, length(tspan));
    % initial time
    t = tspan(1);
    % Propagate the state
    ind = 1;
    while (ind < length(tspan)+1)
        % Changing params
        beta_m = beta_m_0 + omega_e*t;
        xci_m = acos(cos(inc)*cos(gamma_m) + sin(inc)*sin(gamma_m)*cos(Omega - beta_m));
        eta_m = asin(sin(gamma_m)*sin(Omega - beta_m)/sin(xci_m));

        % Magnetic field
        Hb(:, ind) = M/r^3*[cos(thetadot*t - eta_m)*sin(xci_m);
                            cos(xci_m);
                            -2*sin(thetadot*t - eta_m)*sin(xci_m)];

    ind = ind + 1;
    t = t + tstep;
    end
    
end