function [output] = f_n(phi,ellipsoid_name)
%%
% F_N computes the radius of curvature in the prime vertical of a set of 
% points using the geodetic coordinates of the same set of points. The 
% reference ellipsoid should be defined.
%
% HOW: f = f_n(phi,ellipsoid_name)
%
% Input: phi            [n x m] latitude in degrees.
%
%        ellispoid_name         string containing the name of the
%                               reference ellipsoid.
%
% Output: f             [n x m] radius in meters.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: f_w, f_ellipsoid_properties

%% Revision history

%% Remarks

%% Start the algorithm

ell_par = f_ellipsoid_properties(ellipsoid_name); %Define the ellipsoid parameters

output  = ell_par.a./f_w(phi,ellipsoid_name);

end
