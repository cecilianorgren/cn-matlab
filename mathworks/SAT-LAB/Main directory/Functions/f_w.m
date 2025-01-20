function [ output ] = f_w(phi,ellipsoid_name)
%%
% F_W computes the W parameter defined as sqrt(1-(e^2)*(sin(phi)^2) which 
% is used to the computation of the radius of curvature in the prime 
% vertical. The reference ellispoid should be defined.
%
% HOW: f = f_w(phi,ellipsoid_name)
%
% Input: phi            [n x m] latitude in degrees.
%
%        ellispoid_name         string containing the name of the
%                               reference ellipsoid.
%
% Output: f             [n x m] W parameter in meters.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: f_ellipsoid_properties

%% Revision history

%% Remarks

%% Start the algorithm

ell_par = f_ellipsoid_properties(ellipsoid_name); %Define the ellipsoid parameters

output  = sqrt(1-(ell_par.e^2)*(sind(phi).^2));

end
