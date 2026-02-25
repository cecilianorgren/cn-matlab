function out = rotate(in,theta,phi)
% CN.ROTATE Rotates vector in 3D.
%   CN.ROTATE(input,theta,phi) rotates with polar angle theta and 
%   azimuthal angle phi (around 3rd component). Angles given in degrees.
%

% Rotation matrix

R = [cosd(theta) -sind(theta);sind(theta) -cosd(theta)];
switch size(in,2)
    case 3 % no time series
        out = R*in(1:2)';
    case 4 % time series, assuming column one is time
        out = [cosd(theta) -sind(theta);sind(theta) -cosd(theta)]*in(1:2)';
        out = [in(:,1) 
    otherwise
B=B';