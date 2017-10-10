T = 1;
o = 2*pi/T;
nT = 7*T;
nt = 40*nT;
t = linspace(0,nT,nt);
vExB = T/1;
r = 1;

xc = @(t,o,r) r*cos(o*t);
yc = @(t,o,r) r*sin(o*t);

xe = vExB*t + xc(t,o,r);
ye = yc(t,o,r);

xi = vExB*t + xc(t,-o,2*r);
yi = yc(t,-o,2*r)-4;

plot(xe,ye,'linewidth',3,'color',[0.5 0 0.7]); hold on;
plot(xi,yi,'linewidth',3,'color',[1 0.5 0]); hold on;
%arrow([xe(end) ye(end)+0],[xe(end) ye(end)+0.2],'BaseAngle',20)
%arrow([xi(end) yi(end)+0],[xi(end) yi(end)-0.2],'BaseAngle',20)
%quiver(xe(end),ye(end),0,0.1)
axis equal
%axis off
set(gca,'xlim',[-2 10])
hold off;

%%

% Electric field
quiver(-2,-5,0,2,3,'k')
text(-2.7,-0.5,'E','fontsize',20)
% Magnetic field
text(-3,-2,'\circ','fontsize',60)
text(-2.82,-1.4,'.','fontsize',40)
text(-3.5,-1.7,'B','fontsize',20)
%text(-3.2,-2,'o','fontsize',20)

axis off

hold off;
