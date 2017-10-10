function markbeam(varargin)
% mark beam in given or current axis
if nargin == 0
    ax = gca;
else
    ax = varargin{1};
end
irf_pl_mark(ax,toepoch([2007 08 31 10 17 38.158;2007 08 31 10 17 39.689])');