function handles = plot_xo(varargin)

[ax,args,nargs] = axescheck(varargin{:});
plotLabels = 0;
data = args{1};
xyz = args{2};
if nargs>2
  color = args{3};
else 
  color = [0 0 0];%mms_colors('matlab');
end


if isa(data,'TSeries')
  data = data.data;
end
if isa(xyz,'TSeries')
  xyz = xyz.data;
end

if isempty(ax); ax = axes; end
hold(ax,'on')
maxData = max(abs(data));
markerSize = 20;
for iData = 1:size(data,1);
  if data(iData) > 0
    hh = plot(ax,xyz(iData,1),xyz(iData,2),'o','MarkerSize',abs(data(iData))/maxData*markerSize);
  else
    hh = plot(ax,xyz(iData,1),xyz(iData,2),'x','MarkerSize',abs(data(iData))/maxData*markerSize);
  end    
  hh.Color = color;
  hh.LineWidth = 2;
  handles(iData) = hh;
end
hold(ax,'on')



