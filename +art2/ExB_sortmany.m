function [rIndex,zIndex] = ExB_sortmany(x,y,z,rGrid,zGrid,vaztmp)
% Sorts, from contionous x, y and z positions into a r and z predefined
% grid. rGrid and zGrid gives the limits of the slots, so rIndex and zIndex
% ranges from 1 to numel(rGrid)-1.
%
%   [rIndex,zIndex] = ExB_sortmany(x,y,z,rGrid,zGrid);

r = sqrt(x.^2+y.^2);

nr = numel(rGrid-1); % number of r bins
nz = numel(zGrid-1); % number of z bins

rIndex = cell(nr,1);
zIndex = cell(nz,1);

for ii = 1:nz % loop over z bins    
    for jj = 1:nr % loop over r bins
        rlim = rGrid(jj+[0 1]); % edges of r bin
        zlim = zGrid(ii+[0 1]); % edges of z bin
        rind = intersect(find(r>rlim(1)),find(r<rlim(2)));
        zind = intersect(find(z>zlim(1)),find(z<zlim(2)));
        
    end
end        
