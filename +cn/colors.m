function out = colors(n,varargin)
% CN.COLORS    Generates a cell array with different colors.
%       colors = CN.COLORS(n,flag);
%           n - number of colors, currently supporting 6
%               OR, if flag=1 is given, directly pick out that color and
%               return it as a 1x3 RGB array

colors = {[0 0 1],[0 0.8 0],[0 0 0],[0.8 0 0],0.7*[0 1 1],0.8*[1 1 0],0.8*[1 0 1]};
if nargin > 1
    out = colors{n};
else
    out =  colors(1:n);
end