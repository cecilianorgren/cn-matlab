% Plot spacecraft spin plane and E-components
% Sketch, fig 3a in paper

spangle = 30;
rspangle = 0.2;
spanglee = linspace(0,spangle,30);

angle = linspace(0,2*pi,100);
r = 1;
Ereal = 1;
sca = 1.05;

set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

fsize = 16;
% spin plane
plot3(r*cos(angle),r*sin(angle),angle*0,'k'); hold on;
text(-0.2,0.7,0,{'spin','plane'},'fontsize',fsize)

% b-field
bx = cosd(spangle)*[-0.6 1]*1.2;
by = [0 0];
bz = sind(spangle)*[-0.6 1]*1.2;
plot3(bx,by,bz);
text(-0.3,0.2,-0.1,'B','fontsize',fsize)

% angle
plot3(rspangle*cosd(spanglee),spanglee*0,rspangle*sind(spanglee))
text(0.3,0,0.1,'\theta','fontsize',fsize+4)

% e-par
quiver3(0,0,0,Ereal*cosd(spangle),0,Ereal*sind(spangle),sca,'r','linewidth',2)
text(Ereal*cosd(spangle)-0.4,0,Ereal*sind(spangle)-0.0,'E_{||}','fontsize',fsize)

% e-par,sp
quiver(0,0,Ereal*cosd(spangle),0,sca,'g','linewidth',2)
text(Ereal*cosd(spangle)*0.75,-0.3,'E_{||,meas}','fontsize',fsize)

% e-per,sp
quiver(0,0,0,-Ereal*cosd(spangle),sca,'g','linewidth',2)
text(0.1,-0.7,'E_{\perp,meas}','fontsize',fsize)

% projection line
line([Ereal*cosd(spangle) Ereal*cosd(spangle)],[0 0],[0 Ereal*sind(spangle)],'linestyle','--','linewidth',2,'color',[0 0 0])


if 1
    text(-1.25,0,0.92,'a)','fontsize',22)
end




axis equal
axis square
axis off

hold off;