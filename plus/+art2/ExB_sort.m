function [rIndex,zIndex] = ExB_sort(x,y,z,rGrid,zGrid)
% Sorts, from contionous x, y and z positions into a r and z predefined
% grid. rGrid and zGrid gives the limits of the slots, so rIndex and zIndex
% ranges from 1 to numel(rGrid)-1.
%
%   [rIndex,zIndex] = ExB_sort(x,y,z,rGrid,zGrid);

r = sqrt(x.^2+y.^2);
rDifference = rGrid-r;
zDifference = zGrid-z;

rIndex = find(abs(rDifference)==min(abs(rDifference)));
zIndex = find(abs(zDifference)==min(abs(zDifference)));