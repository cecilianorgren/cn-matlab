function [value,isterminal,direction] = events(t,y,L)
% Locate the time when the particle exits the box by top or bottom and stop
% integration.

value = 0.5*L*1e3 - abs(y(3))*0.999; % detect z>L/2 and z<L/2 (value = 0)
isterminal = 1;   % Stop the integration
direction = 0;   % Both directions