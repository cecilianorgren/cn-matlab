classdef f_ellipsoid
%%
% F_ELLIPSOID defines a new class with the fundamental parameters of an
% ellipsoid.
%
% HOW: f = f_ellipsoid
%
% Output: f            [1 x 1] element defined as ellipsoid.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

    properties
        
        epsg_code   %Code of EPSG/OGP Geodetic Parameter Database
        a           %Semi-Major axis
        b           %Semi-Minor axis
        f           %Flattening
        f_inv       %Inverse flattening
        e           %First eccentricity
        e_prime     %Second eccentricity
        gm          %Gravitational constant
        omega       %Angular velocity
        
    end
    
end
