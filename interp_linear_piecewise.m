function new_data = interp_linear_piecewise(old_data,old_time,new_time)
% INTERP_LINEAR_PIECEWISE Connecting the dots.
%   new_data = interp_linear_piecewise(old_data,old_time,new_time);
%   new_data 

size_old_data = size(old_data);
size_new_data = size_old_data;
size_new_data(size_new_data>1) = numel(new_time);

old_data = torow(old_data);
old_time = torow(old_time);
new_time = torow(new_time);

[~,k] = histc(new_time,old_time); % k is the bin
n = length(old_time);
k(k == n) = n - 1;
t = (new_time - old_time(k))./(old_time(k+1) - old_time(k));
new_data = (1-t).*old_data(k) + t.*old_data(k+1);

new_data = reshape(new_data,size_new_data);