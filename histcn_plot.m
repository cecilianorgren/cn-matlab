function out = histcn_plot(varargin)

[N edges mid loc] = histcn(varargin{:});
pcolor(mid{1:2},log10(N(:,:,ceil(end/2)))')
shading('flat')
colorbar
