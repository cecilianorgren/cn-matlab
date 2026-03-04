function [i,j,k] = find_peaks_iteratively(F,N,w)
% Find peaks iteratively in 3D matrix.
%   1. Find max.
%   2. Remove 

npeaks = 2; % assume max 2 cold distinct peaks
i = zeros(npeaks,1);
j = zeros(npeaks,1);
k = zeros(npeaks,1);
w = 3;
for ipeak = 1:npeaks
  [maxval, idx] = max(Fiter(:));
  [i(ipeak),j(ipeak),k(ipeak)] = ind2sub(size(Fiter), idx);
  Fiter(i(ipeak)+(-w:1:w),j(ipeak)+(-w:1:w),k(ipeak)+(-w:1:w)) = 0;
end
