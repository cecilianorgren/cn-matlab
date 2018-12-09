function [ind1 ind2] = ind12(n,nint)
% divides 
n_per_interval = floor(n/nint);
ind1 = 1:n_per_interval:(n-1);
ind2 = [ind1(2:end)-1 n];