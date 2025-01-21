function [output_x,output_y,output_z] = f_geo2cart(latitude,longitude,h,ellipsoid_name)
%%
% F_GEO2CART computes the cartesian coordinates of a set of points using 
% the geodetic coordinates of the same set of points. The reference 
% ellispoid should be defined.
%
% HOW: [x,y,z] = f_geo2cart(phi,lamda,h,ellipsoid_name)
%
% Input: latitude       [n x m] latitude in degrees.
%
%        longitude      [n x m] longitude in degrees.
%
%        h              [n x m] ellipsoidal height in meters.
%
%        ellispoid_name         string containing the name of the
%                               reference ellipsoid.
%
% Output: x             [n x m] x coordinate in meters.
%
%         y             [n x m] y coordinate in meters.
%
%         z             [n x m] z coordinate in meters.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: f_n, f_ellipsoid_properties

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 4
    error('Wrong number of input arguments')
end

if max(size(latitude) ~= size(longitude)) == 1
    error('phi and lamda do not have the same dimensions')
end

ell_par  = f_ellipsoid_properties(ellipsoid_name); %Define the ellipsoid parameters

n        = f_n(latitude,ellipsoid_name);

output_x = (n+h).*cosd(latitude).*cosd(longitude); %X-coordinate
output_y = (n+h).*cosd(latitude).*sind(longitude); %Y-coordinate
output_z = ((1-ell_par.e^2).*n+h).*sind(latitude); %Z-coordinate

end
