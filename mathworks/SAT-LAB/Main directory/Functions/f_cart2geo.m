function [output_latitude,output_longitude,output_h] = f_cart2geo(x,y,z,ellipsoid_name)
%%
% F_CART2GEO computes the geodetic coordinates of a set of points using 
% the cartesian coordinates of the same set of points. The reference 
% ellispoid should be defined.
%
% HOW: [latitude,longitude,h] = f_cart2geo(x,y,z,ellipsoid_name)
%
% Input: x              [n x m] x coordinate in meters.
%
%        y              [n x m] y coordinate in meters.
%
%        z              [n x m] z coordinate in meters.
%
%        ellispoid_name         string containing the name of the
%                               reference ellipsoid.
%
% Output: latitude      [n x m] latitude in degrees.
%
%         longitude     [n x m] longitude in degrees.
%
%         h             [n x m] ellipsoidal height in meters.
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

if max(size(x) ~= size(y) ) == 1 || max(size(y) ~= size(z) ) == 1
    error('Input coordinates do not have the same dimensions')
end

%Define the ellipsoid parameters
ell_par             = f_ellipsoid_properties(ellipsoid_name); 

%Initialize phi
output_latitude     = atand(z.*(1+(ell_par.e_prime^2))./sqrt(x.^2+y.^2));                                                    %Degrees

%Iterations for phi
for i = 1:10
    
    output_latitude = atand((z+(ell_par.e^2).*f_n(output_latitude,ellipsoid_name).*sind(output_latitude))./sqrt(x.^2+y.^2)); %Degrees
    
end

output_longitude    = f_atan2d(y,x);                                                                                         %Degrees
output_h            = (z./sind(output_latitude))-(1-ell_par.e^2)*f_n(output_latitude,ellipsoid_name);                        %Meters

end
