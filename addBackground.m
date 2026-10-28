function h = addBackground(ax,fun,val,varargin)


for iax = 1:numel(ax)
  x = ax(iax).XLim;
  y = ax(iax).YLim;
  x = linspace(x(1),x(2),100);
  y = linspace(y(1),y(2),101);
  [X,Y] = ndgrid(x,y);
  C = fun(X,Y);
  hold(ax(iax),'on')
  contour(ax(iax),X,Y,C,val,varargin{:})
  %pcolor(ax(iax),X,Y,C)
  hold(ax(iax),'off')
end