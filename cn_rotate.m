function out = cn_rotate(in,theta,phi)
% Rotates vector with polar angle theta and azimuthal angle phi
R = [cosd(theta) -sind(theta);sind(theta) -cosd(theta)];
switch size(in,2)
    case 3 % no time series
        out = R*in(1:2)';
    case 4 % time series, assuming column one is time
        out = [cosd(theta) -sind(theta);sind(theta) -cosd(theta)]*in(1:2)';
        out = [in(:,1) 
else
B=B';