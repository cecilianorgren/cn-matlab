%% geometric lines background

nPoints = 10;
xMax = 1;
yMax = 1;
zMax = 1;
x = xMax*(1 + 2*randn(nPoints,1));
y = yMax*(1 + 2*randn(nPoints,1));
z = zMax*(1 + 2*randn(nPoints,1));
%x = xMax*randn(nPoints,1);
%y = yMax*randn(nPoints,1);
%z = zMax*randn(nPoints,1);




cMap = colormap('jet');
cPoints = [interp1(1:size(cMap,1),cMap(:,1),linspace(1,size(cMap,1),nPoints))',...
           interp1(1:size(cMap,1),cMap(:,2),linspace(1,size(cMap,1),nPoints))',...
           interp1(1:size(cMap,1),cMap(:,3),linspace(1,size(cMap,1),nPoints))'];
hca = subplot(1,1,1);
scatter3(hca,x,y,z,10*ones(nPoints,1),1:nPoints)
hold(hca,'on')
for iPoint1 = 1:nPoints
  for iPoint2 = 1:nPoints
    cPoint = mean([cPoints(iPoint1,:);cPoints(iPoint2,:)],1);
    %plot(hca,[x(iPointX) x(iPointY)],[y(iPointX) y(iPointY)],'color',cPoint)
    plot3(hca,[x(iPoint1) x(iPoint2)],[y(iPoint1) y(iPoint2)],[z(iPoint1) z(iPoint2)],'color',cPoint)
  end
end
hold(hca,'off')
axis(hca,'off')
set(gcf,'color','white');

%% geometric lines background

nPoints = 10;
xMax = 1;
yMax = 1;
zMax = 1;
x = xMax*(1 + 2*randn(nPoints,1));
y = yMax*(1 + 2*randn(nPoints,1));
z = zMax*(1 + 2*randn(nPoints,1));
x = xMax*randn(nPoints,1);
y = yMax*randn(nPoints,1);
z = zMax*randn(nPoints,1);

cMap = colormap('jet');
cPoints = [interp1(1:size(cMap,1),cMap(:,1),linspace(1,size(cMap,1),nPoints))',...
           interp1(1:size(cMap,1),cMap(:,2),linspace(1,size(cMap,1),nPoints))',...
           interp1(1:size(cMap,1),cMap(:,3),linspace(1,size(cMap,1),nPoints))'];
hca = subplot(1,1,1);
scatter3(hca,x,y,z,10*ones(nPoints,1),1:nPoints)
hold(hca,'on')
for iPoint1 = 1:nPoints
  for iPoint2 = 1:nPoints
    cPoint = mean([cPoints(iPoint1,:);cPoints(iPoint2,:)],1);
    %plot(hca,[x(iPointX) x(iPointY)],[y(iPointX) y(iPointY)],'color',cPoint)
    plot3(hca,[x(iPoint1) x(iPoint2)],[y(iPoint1) y(iPoint2)],[z(iPoint1) z(iPoint2)],'color',cPoint)
  end
end
hold(hca,'off')

