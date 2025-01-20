function R_z = f_rotz(theta)
%%
% F_ROTZ constructs the 3D rotation matrix for rotations around z-axis.
%
% HOW: R_z = f_rotz(theta)
%
% Input: theta          [1 x 1] rotation angle in rad.
%
% Output: R_z           [3 x 3] rotation matrix.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 23/01/2017
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 1
    error('Wrong number of input arguments')
end

if max(size(theta)) ~= 1
    error('Select only one rotation angle')
end

R_z = [ cos(theta), sin(theta), 0 ;
       -sin(theta), cos(theta), 0 ;
        0         , 0         , 1];

end
