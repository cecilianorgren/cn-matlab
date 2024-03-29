nlines = 5;      

colors = [...
    0    0.4470    0.7410 ;...
    0.8500    0.3250    0.0980 ;...
    0.9290    0.6940    0.1250 ;...
    0.4940    0.1840    0.5560 ;...
    0.4660    0.6740    0.1880 ;...
    0.3010    0.7450    0.9330 ;...
    0.6350    0.0780    0.1840];
  
r = @(v,phi,phi0,w,r0) v*(phi-phi0)/w+r0;

r0 = 4;
r0dot = 0.5;
v = 0.1;
phi0 = -60;
phiend = 00;
w = 0.1;

phi = linspace(phi0,phiend,nlines);
philine = linspace(phi0,phiend,100);

%close all
hca = subplot(1,1,1);
% Sun
hp = patch(hca,r0*cosd(linspace(0,360,100)),r0*sind(linspace(0,360,100)),colors(3,:));
hold(hca,'on')

% Spiral
rspiral = r(v,philine,phi0,w,r0);
xxspiral = rspiral.*cosd(philine);
yyspiral = rspiral.*sind(philine);
plot(hca,xxspiral,yyspiral,'linewidth',1,'color',[1 0 0])

% Radial lines
linemarks = nlines:-1:1;
for iphi = 1:nlines  
  phitmp = phi(iphi);
  rtmp = r(v,phitmp,phi0,w,r0);
  xx = [r0 rtmp]*cosd(phitmp);
  yy = [r0 rtmp]*sind(phitmp);
  plot(hca,xx,yy,'linewidth',1,'color',[0 0 0])
  patch(hca,xx(2)+r0dot*cosd(linspace(0,360,100)),...
           yy(2)+r0dot*sind(linspace(0,360,100)),colors(1,:))
  %text(hca,xx(2),yy(2),sprintf('t_%g',linemarks(iphi)),'horizontalalignment','left')
end
hold(hca,'off')

axis(hca,'equal')
hca.YDir = 'reverse';
