function h = addBackground(ax,fun,val,varargin)
% ADDBACKGROUND draws equicontourlines for the pased function
% 
%   Examples:
%   Straight lines: y = k*x -> y./x = k
%   addBackground(ax,@(x,y)y./x,k);
%   
%   Circular shells: x^2 + y^2 = r^2 -> sqrt(x.^2+y.^2) = r
%   addBackground(ax,@(x,y)sqrt(x.^2+y.^2),r);

for iax = 1:numel(ax)
  clim = ax(iax).CLim;
  x = ax(iax).XLim;
  y = ax(iax).YLim;
  x = linspace(x(1),x(2),1000);
  y = linspace(y(1),y(2),1001);
  [X,Y] = ndgrid(x,y);
  C = fun(X,Y);
  hold(ax(iax),'on')
  if numel(val) == 1; val = [val val]; end
  contour(ax(iax),X,Y,C,val,varargin{:},'color','k')
  %pcolor(ax(iax),X,Y,C)
  hold(ax(iax),'off')
  ax(iax).CLim = clim;
end
