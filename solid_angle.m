function out = solid_angle(polar_angle,azimuthal_angle)
% Calculates solid angle int sin(theta) dtheta dphi
%   out = solid_angle(polar_angle,azimuthal_angle)

out = diff(azimuthal_angle,1,2).*(cosd(polar_angle(:,1))-cosd(polar_angle(:,2)));
out = out*pi/180;



