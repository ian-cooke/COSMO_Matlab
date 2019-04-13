% tspan is time span for integration
% y0 is initial state vector, 9x1
% I is the intertia matrix, 3x3

function [t, y, omegaDot] = integrateMRPs(tspan, y0, I)

t = tspan(1);
y = y0';

% opts = odeset('Events', @mrpODEevents, 'OutputFcn', @mrpODEOutputFcn);
opts = odeset('Events', @mrpODEevents);
te = tspan(1);
ie = 1;


c = 1;
while ie == 1
    %     [tTemp,yTemp,te,ye,ie] = ode45(@(t,y) mrpODE(t,y,diag(I)), [te(end) tspan(2)],y0, opts);
    [tTemp,yTemp,te,ye,ie] = ode45(@(t,y) mrpODE(t,y,I), [te(end) tspan(2)],y0, opts);
    if ie == 1
        y0 = ye(end,:)';
        sig = y0(1:3);
        sig = -sig/(sig'*sig);
        y0(1:3) = sig;
    end
    c = c+1;
    
    y = [y;yTemp];
    t = [t;tTemp];
end

omegaDot = zeros(length(y),3);
omegaDot(1,:) = [0 0 0];

for i = 2:length(y)
    omegaDot(i,:) = (y(i,4:6) - y(i-1,4:6))/(t(i) - t(i-1));
    for j = 1:3
        if isnan(omegaDot(i,j))
            omegaDot(i,j) = 0;
        end
    end
end



