function [dydt wDotT] = mrpODE(t,y,I)

% persistent wDot
% global omegaDot

% if nargout > 1
%     dydt = [];
%     wDotT = wDot;
%     return
% end

% theta = y(1:3);
sigma = y(1:3);
omega = y(4:6); % thetaDot

sigmaTilde = makeSkew(sigma);
dsigma_dt = 0.25 * [(1-(sigma'*sigma))*eye(3) + 2*sigmaTilde + 2*(sigma*sigma')]*omega;

omegaDot = inv(I)*(-makeSkew(omega)*I*omega);

% wDot = [wDot; omegaDot' t];

dydt(1:3) = dsigma_dt;
dydt(4:6) = omegaDot;
dydt = dydt';


