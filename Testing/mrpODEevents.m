function [position,isterminal,direction] = mrpODEevents(t,y)
position = norm(y(1:3)) - 1 - 1e-5; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction