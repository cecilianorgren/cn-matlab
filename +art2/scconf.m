close

cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields
if ~exist('ind','var'); ind = 5; end
tint = tint{highQuality(ind)};

% dp = R4-R3, -> R4 = R3 + dp
dx=cn.mean(irf_tlim(dp1,tint),1); % SP
dy=cn.mean(irf_tlim(dp2,tint),1); % BxSP
dz=cn.mean(irf_tlim(dp3,tint),1); % B

R3 = [0 0 0];
R4 = R3 + [dx dy dz];

xx = [R3(1) R4(1) R4(1) R4(1)];
yy = [R3(2) R3(2) R4(2) R4(2)];
zz = [R3(3) R3(3) R3(3) R4(3)];

line(xx,yy,zz,'color',[0 0 0]); hold on;
h = plot3(R3(1),R3(3),R3(3),'go',R4(1),R4(2),R4(3),'bo');
set(h,'markersize',14,'linewidth',2,'marker','square')%,'ydata',y*1.2*[-1 1],'zdata',z*1.2*[-1 1])
text(R3(1),R3(3),R3(3),'      C3')
text(R4(1),R4(2),R4(3),'      C4')

text(0.5*R4(1),0*R4(2),0*R4(3),[num2str(dx,'%.0f') ' km'])
text(0.5*R4(1),0.5*R4(2),0*R4(3),[num2str(dy,'%.0f') ' km'])
text(0.5*R4(1),1*R4(2),0.5*R4(3),[num2str(dz,'%.0f') ' km'])

view([1 -1 1])
axis equal
axis off

% add spin plane tilt



hold off


