function hq = plot_quivers_tmp(varargin)
% plot_quivers
%   plot_quivers(ax,data,xyz)
%   plot_quivers(data,xyz)

[ax,args,nargs] = axescheck(varargin{:});
plotLabels = 0;
data = args{1};
xyz = args{2};
S = 0;
hq = [];

if nargs>2
  S = args{3}; 
end
if nargs>3
  color = args{4};
else 
  color = [0 0 0];
end
  
if nargs>4
  labels = args{4};
  plotLabels = 1;
end


if isa(data,'TSeries')
  data = data.data;
end
if isa(xyz,'TSeries')
  xyz = xyz.data;
end

if size(data,2) == 2
%   if isa(data,'TSeries')
%     qx = data.x.data;
%     qy = data.y.data;
%   else
    qx = data(:,1);
    qy = data(:,2);
%   end
%   if isa(xyz,'TSeries')
%     x = xyz.x.data;
%     y = xyz.y.data;
%   else
    x = xyz(:,1);
    y = xyz(:,2);
%   end    
    if isempty(ax)
      hq = quiver(x,y,qx,qy,S,'color',color,'linewidth',1.5);
    else
      hq = quiver(ax,x,y,qx,qy,S,'color',color,'linewidth',2);
    end
elseif size(data,2) == 3  
%   if isa(data,'TSeries')
%     qx = data.x.data;
%     qy = data.y.data;
%     qz = data.z.data;
%   else
    qx = data(:,1);
    qy = data(:,2);
    qz = data(:,3);
%   end
%   if isa(xyz,'TSeries')
%     x = xyz.x.data;
%     y = xyz.y.data;
%     z = xyz.z.data;
%   else
    x = xyz(:,1);
    y = xyz(:,2);
    z = xyz(:,3);
%   end
  % permute data order
%   newx = z;
%   newy = -y;
%   newz = x;
%   x = newx;
%   y = newy;
%   z = newz;
  if plotLabels
    for iText = 1:numel(labels) 
      text(double(x(iText)+qx(iText)),double(y(iText)+qy(iText)),double(z(iText)+qz(iText)),labels{iText})
    end
  end
  if isempty(ax)
    quiver3(x,y,z,qx,qy,qz,'color',color)
  else
    quiver3(ax,x,y,z,qx,qy,qz,'color',color)
  end
end






%set(ax,'xdir','reverse');
%set(ax,'ydir','reverse');

