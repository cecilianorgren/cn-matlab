function plot_quivers(varargin)
% plot_quivers
%   plot_quivers(ax,data,xyz)
%   plot_quivers(data,xyz)

[ax,args,nargs] = axescheck(varargin{:});
plotLabels = 0;
data = args{1};
xyz = args{2};
if nargs>2
  color = args{3};
else 
  color = [0 0 0];
end
if nargs>3
  labels = args{4};
  plotLabels = 1;
end
  
S = 0.5;


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
      quiver(x,y,qx,qy,S,'color',color)
    else
      quiver(ax,x,y,qx,qy,S,'color',color)
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
    quiver3(x,y,z,qx,qy,qz,S,'color',color)
  else
    quiver3(ax,x,y,z,qx,qy,qz,S,'color',color)
  end
end






%set(ax,'xdir','reverse');
%set(ax,'ydir','reverse');

