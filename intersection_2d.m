function varargout = intersection_2d(X,Y,mapA,mapB,levA,levB,varargin)
% INTERSECTION_2D Calculate intersection of 2D equicontour levels 
% 

% Default values
doTSeries = 0;
doPlot = 0;
doMesh = 0;

% Collect input
[~,args,nargs] = axescheck(varargin{:});
nargin = numel(varargin);

have_options = nargs > 1;
while have_options
  switch(lower(args{1}))
    case 'plot' % velocity interval to consider
      l = 2;
      doPlot = args{2};      
    case 'mesh' % velocity interval to consider
      l = 2;
      doMesh = args{2};
    otherwise
      fprintf('Unknown input argument: %s.\n',args{1})
      l = 1;
  end
  args = args((l+1):end);
  if isempty(args), break, end
end

%% Calculate intersection between lines
if isa(levA,'TSeries')
  levA_orig = levA;
  levA = double(levA.data);
  levB_orig = levB;
  levB = double(levB.data);
  doTSeries = 1;
end
if doMesh
  [levA,levB] = ndgrid(levA,levB);
end
nLevelsX = size(levA,1);
nLevelsY = size(levA,2);

% Are these changing order in some way? Perhaps do within loop?
%Cx = contourcs(X(:,1),Y(1,:),mapA',levA);
%Cy = contourcs(X(:,1),Y(1,:),mapB',levB);


% Use polyxpoly to find intersections

x_intersect = nan(nLevelsX,nLevelsY);
y_intersect = nan(nLevelsX,nLevelsY);

for ix = 1:nLevelsX
  for iy = 1:nLevelsY
    % Check if level value is within map range
    if levA(ix,iy)>=min(mapA(:)) && levA(ix,iy)<=max(mapA(:)) && levB(ix,iy)>=min(mapB(:)) && levB(ix,iy)<=max(mapB(:))
      CA = contourcs(X(:,1),Y(1,:),mapA',levA(ix,iy)*[1 1]);
      CB = contourcs(X(:,1),Y(1,:),mapB',levB(ix,iy)*[1 1]);       
    else
      continue
    end

    [xint_tmp, yint_tmp] = polyxpoly(CA(1).X, CA(1).Y, CB(1).X, CB(1).Y);               
    
    try
    if not(isempty(xint_tmp))
      if numel(xint_tmp) == 1
        x_intersect(ix,iy) = xint_tmp;
        y_intersect(ix,iy) = yint_tmp;
      else
        x_intersect(ix,iy) = mean(xint_tmp);
        y_intersect(ix,iy) = mean(yint_tmp);
        disp(sprintf('Numerous crossings (%g) at: [%g,%g]',numel(xint_tmp),levA(ix,iy),levB(ix,iy)))
      end
    end
    catch
      1;
    end
  end
end
try
if doPlot
  %%
  hca = subplot(1,2,1); 
  [c,hcc] = contour(hca,X,Y,mapA,unique(levA(isfinite(levA))),'k');
  %clabel(c,hcc)
  hold(hca,'on')
  contour(hca,X,Y,mapB,unique(levB(isfinite(levB))),'r')
  plot(hca,x_intersect,y_intersect,'bo')
  hold(hca,'off')
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  
  hca = subplot(1,2,2); 
  [c,hcc] = contour(hca,X,Y,mapA,0.1:0.1:1,'k');
  clabel(c,hcc);
  hold(hca,'on')
  [c,hcc] = contour(hca,X,Y,mapB,0.2:0.2:2,'r');
  clabel(c,hcc);
  plot(hca,x_intersect,y_intersect,'bo')
  hold(hca,'off')
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
end
catch
  1;
end
if doTSeries
  x_intersect = irf.ts_scalar(levA_orig.time,x_intersect);
  y_intersect = irf.ts_scalar(levA_orig.time,y_intersect);
end
varargout{1} = x_intersect;
varargout{2} = y_intersect;
