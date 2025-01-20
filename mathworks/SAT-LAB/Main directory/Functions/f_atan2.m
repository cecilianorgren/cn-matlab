function [output] = f_atan2(dx,dy)
%%
% F_ATAN2 computes the four-quadrant inverse tangent (see also MATLAB atan2
% function). The result is given in rad.
%
% HOW: f = f_atan2(dx,dy)
%
% Input: dx            [n x m] numerator.
%
%        dy            [n x m] denominator.
%
% Output: f            [n x m] four-quadrant inverse tangent in rad.
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
    
    dy                = dx(:,2);
    dx                = dx(:,1);
    
end

if max(size(dx) ~= size(dy)) == 1
    error('dx and dy do not have the same dimensions')
end

output                = NaN(size(dx));

g                     = atan(abs(dx./dy)); %Two-quadrant inverse tangent

output(dx>0  & dy>0)  = g(dx>0 & dy>0);
output(dx>0  & dy==0) = pi/2;
output(dx>0  & dy<0)  = pi-g(dx>0 & dy<0);
output(dx==0 & dy>0)  = 0;
output(dx==0 & dy<0)  = pi;
output(dx<0  & dy>0)  = 2*pi-g(dx<0 & dy>0);
output(dx<0  & dy==0) = 3*pi/2;
output(dx<0  & dy<0)  = pi+g(dx<0 & dy<0);

end
