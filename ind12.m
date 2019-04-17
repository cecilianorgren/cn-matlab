function [ind1 ind2] = ind12(n,nint)
% divides 
rest = mod(n,nint);
n_per_interval = floor(n/nint);
ind1 = 1:n_per_interval:(n-rest);
ind2 = [ind1(2:end)-rest n];