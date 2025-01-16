function out = histcn_plot(varargin)

%[ax,args,nargs] = axescheck(varargin);
%args = args{1};
args = varargin;
[N edges mid loc] = histcn(args{:});

%if isempty(ax)
ax = gca;
%end

pcolor(ax,mid{1:2},log10(N(:,:,ceil(end/2)))')
shading(ax,'flat')
colorbar(ax)
