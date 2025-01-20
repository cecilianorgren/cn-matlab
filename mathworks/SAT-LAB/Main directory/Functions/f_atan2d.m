function [output] = f_atan2d(dx,dy)
%%
% F_ATAN2D computes the four-quadrant inverse tangent (see also MATLAB
% atan2d function). The result is given in degrees.
%
% HOW: f = f_atan2d(dx,dy)
%
% Input: dx            [n x m] numerator.
%
%        dy            [n x m] denominator.
%
% Output: f            [n x m] four-quadrant inverse tangent in degrees.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin < 1 || nargin > 2
    error('Wrong number of input arguments')
elseif nargin == 1
    
    dy = dx(:,2);
    dx = dx(:,1);
    
end

if max(size(dx) ~= size(dy)) == 1
    error('dx and dy do not have the same dimensions')
end

output = f_atan2(dx,dy)*180/pi;

end
