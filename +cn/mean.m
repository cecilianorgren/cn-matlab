function out = mean(in,varargin)
% CN.MEAN   Takes the mean of a time series.
%       out = CN.MEAN(in,flag) assumes in is a time series where the first
%       column is the time and the rest are data. If in has dimensions NxM,
%       then out has dimensions 1xM, where the first column is the mean 
%       time. If flag=1 is given, then out is simply the value, without
%       time indication, 1x(M-1).

if nargin>1 && varargin{1}==1
    out = mean(in(:,2:end,1));
else
    out = mean(in);
end

