%% 2d two structures directly
fig=figure;
set(gcf,'PaperPositionMode','auto');
    
if 1 % Converging E 
[X1E,Y1E] = meshgrid(-2:.4:2); Z1E = 2.*exp(-X1E.^2 - Y1E.^2);
[DX1,DY1] = gradient(Z1E,.4,.4);

[X1phi,Y1phi] = meshgrid(-2:.1:2); Z1phi = 1.*exp(-X1phi.^2 - Y1phi.^2);

colormap hsv
cont1=contour(X1phi,Y1phi,Z1phi,3); % phi


hold on
q1=quiver(X1E,Y1E,DX1,DY1,0.5,'color',[0 0.75 0.75]); set(q1,'linewidth',1) % E

icol=[0.2 0.7 0.2];
r=0.5;
ang=-pi/4:0.001:pi/4-0.1;
hold on;
re=0.2;
vex=cos(ang)+re*cos(ang*20);
vey=sin(ang)+re*sin(ang*20);
plot3(vex,vey,ang*0,'linewidth',1.1,'color',[0 0 0]) % e-

iang=0:0.01:2*pi;
nang=size(iang,2);
or=[0 0 -1.5];
ix=or(1)+cos(iang);
iy=or(2)+sin(iang);
plot3(ix,iy,ix*0+or(3),'linewidth',1.1,'color',[0.2 0.7 0.2]) % i

b1col=[1 0 0];
b0=quiver3(0,0,0,0,0,1,1,'linewidth',1.1,'color',[0.6 0 0]) % b0
b1=quiver3(or(1),or(2),or(3),0,0,-0.5,1,'linewidth',1.1,'color',[1 0 0]) % b1
%b1arrow=arrow([or(1) or(2) or(3)], [or(1) or(2) or(3)-0.1],...
%    'Length',1,'Width',10,...
%    'FaceColor',b1col,'EdgeColor',b1col)

vearrow=arrow([vex(end) vey(end)], [vex(end)-0.1 vey(end)+0.1],...
    'Length',30)
iarrow1=arrow([ix(end) iy(end)-0.1 or(3)], [ix(end) iy(end)+0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol)
%set('color',[1 0 1])
iarrow2=arrow([ix(fix(nang/2)) iy(fix(nang/2))+0.01 or(3)],...
    [ix(fix(nang/2)) iy(fix(nang/2))-0.1 or(3)],...
    'Length',20,'width',10,...
    'FaceColor',icol,'EdgeColor',icol)
end

if 1 % Diverging E 
[X,Y] = meshgrid(2:.4:6,-2:.4:2);
Z = -1.*exp(-(X-4).^2 - Y.^2);
[DX,DY] = gradient(Z,.4,.4);

[X2,Y2] = meshgrid(2:.1:6,-2:.1:2);
Z2 = -1.*exp(-(X2-4).^2 - Y2.^2);

colormap hsv
cont2=contour(X2,Y2,Z2,4); % phi
hold on
q2=quiver(X,Y,DX,DY,0.5); set(q1,'linewidth',1) % E

r=0.5;
ang=-pi/4:0.001:pi/4;
hold on;
re=0.2;
vexb=4+cos(-ang)+re*cos(ang*20);
veyb=+sin(-ang)+re*sin(ang*20);
plot3(vexb,veyb,ang*0,'linewidth',1.1,'color',[0 0 0]) % e-

icol=[0.2 0.7 0.2];
iang=0:0.01:2*pi;
nang=size(iang,2);
or=[4 0 -1.5];
ix=or(1)+cos(iang);
iy=or(2)+sin(iang);
plot3(ix,iy,ix*0+or(3),'linewidth',1.1,'color',icol) % i

b0b=quiver3(4,0,0,0,0,1,1,'linewidth',1.1,'color',[0.7 0 0]) % b0
b1b=quiver3(or(1),or(2),or(3),0,0,0.5,1,'linewidth',1.1,'color',[1 0 0]) % b1


vearrow=arrow([vexb(end) veyb(end)], [vexb(end)-0.1 veyb(end)-0.1],...
    'Length',30)
iarrow1=arrow([ix(end) iy(end)+0.1 or(3)], [ix(end) iy(end)-0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol)
%set('color',[1 0 1])
iarrow2=arrow([ix(fix(nang/2)) iy(fix(nang/2))-0.1 or(3)],...
    [ix(fix(nang/2)) iy(fix(nang/2))+0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol)
end
%linkaxes(q1,q2)
%colorbar
hold off
axis equal
legend('\phi','E','e^-','I','B_0','B_1')
%etxt=text(1.2,0,0.5,'e','color',[0 0 0],'fontsize',20);


set(gca,'ylim',[-3 3],'xlim',[-2 6],'zlim',[-2 1],...
    'ycolor','white','xcolor',[1 1 1],'zcolor',[1 1 1]);
set(gcf,'color',[1 1 1])
set(gca,'CameraPosition',[2.8044  -46.0842   24.0071])
eval(['print -depsc2 db_sketchie.eps']);

%% 2d (the one I use)
fig=figure;
set(fig,'defaulttextfontsize',10)
set(gcf,'PaperPositionMode','auto');
if 1 % Converging E 
[X,Y] = meshgrid(-2:.4:2);
Z = 1.*exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,.4,.4);

[X2,Y2] = meshgrid(-2:.1:2);
Z2 = 1.*exp(-X2.^2 - Y2.^2);

colormap hsv
cont=contour(X2,Y2,Z2,3); % phi


hold on
q1=quiver(X,Y,DX,DY,0.5,'color',[0 0.75 0.75]); set(q1,'linewidth',1) % E

icol=[0.2 0.7 0.2];
r=0.5;
ang=-pi/4:0.001:pi/4-0.1;
hold on;
re=0.2;
vex=cos(ang)+re*cos(ang*20);
vey=sin(ang)+re*sin(ang*20);
plot3(vex,vey,ang*0,'linewidth',1.1,'color',[0 0 0]) % e-

iang=0:0.01:2*pi;
nang=size(iang,2);
or=[0 0 -1.5];
ix=or(1)+cos(iang);
iy=or(2)+sin(iang);
plot3(ix,iy,ix*0+or(3),'linewidth',1.1,'color',[0.2 0.7 0.2]) % i

b1col=[1 0 0];
b0=quiver3(0,0,0,0,0,1,1,'linewidth',1.1,'color',[0.6 0 0]); % b0
b1=quiver3(or(1),or(2),or(3),0,0,-0.4,1,'linewidth',1.1,'color',[1 0 0]); % b1
%b1arrow=arrow([or(1) or(2) or(3)], [or(1) or(2) or(3)-0.1],...
%    'Length',1,'Width',10,...
%    'FaceColor',b1col,'EdgeColor',b1col)

vearrow=arrow([vex(end) vey(end)], [vex(end)-0.1 vey(end)+0.1],...
    'Length',30);
iarrow1=arrow([ix(end) iy(end)+0.1 or(3)], [ix(end) iy(end)-0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol);
%set('color',[1 0 1])
iarrow2=arrow([ix(fix(nang/2)) iy(fix(nang/2))-0.1 or(3)],...
    [ix(fix(nang/2)) iy(fix(nang/2))+0.1 or(3)],...
    'Length',20,'width',10,...
    'FaceColor',icol,'EdgeColor',icol);
end

if 1 % Diverging E 
[X,Y] = meshgrid(2:.4:6,-2:.4:2);
Z = -1.*exp(-(X-4).^2 - Y.^2);
[DX,DY] = gradient(Z,.4,.4);

[X2,Y2] = meshgrid(2:.1:6,-2:.1:2);
Z2 = -1.*exp(-(X2-4).^2 - Y2.^2);

colormap hsv
cont2=contour(X2,Y2,Z2,3); % phi
hold on
q2=quiver(X,Y,DX,DY,0.5); set(q1,'linewidth',1) % E

r=0.5;
ang=-pi/4:0.001:pi/4;
hold on;
re=0.2;
vexb=4+cos(-ang)+re*cos(ang*20);
veyb=+sin(-ang)+re*sin(ang*20);
plot3(vexb,veyb,ang*0,'linewidth',1.1,'color',[0 0 0]) % e-

icol=[0.2 0.7 0.2];
iang=0:0.01:2*pi;
nang=size(iang,2);
or=[4 0 -1.5];
ix=or(1)+cos(iang);
iy=or(2)+sin(iang);
plot3(ix,iy,ix*0+or(3),'linewidth',1.1,'color',icol) % i

b0b=quiver3(4,0,0,0,0,1,1,'linewidth',1.1,'color',[0.7 0 0]); % b0
b1b=quiver3(or(1),or(2),or(3),0,0,0.4,1,'linewidth',1.1,'color',[1 0 0]); % b1


b1col=[1 0 0];
b0col=[0.7 0 0];

vearrow=arrow([vexb(end) veyb(end)], [vexb(end)-0.1 veyb(end)-0.1],...
    'Length',30);
iarrow1=arrow([ix(end) iy(end)-0.1 or(3)], [ix(end) iy(end)+0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol);
%set('color',[1 0 1])
iarrow2=arrow([ix(fix(nang/2)) iy(fix(nang/2))+0.1 or(3)],...
    [ix(fix(nang/2)) iy(fix(nang/2))-0.1 or(3)],...
    'Width',10,'length',20,...
    'FaceColor',icol,'EdgeColor',icol);
 db1a=arrow([4	-0.1 or(3)+0.4], [4 0.1 or(3)+0.4],...
     'Width',10,'length',20,...
     'FaceColor',b1col,'EdgeColor',b1col);
 db1b=arrow([0	0.1 or(3)-0.4], [0 -0.1 or(3)-0.4],...
     'Width',10,'length',20,...
     'FaceColor',b1col,'EdgeColor',b1col);
 db0a=arrow([4	-0.1 1], [4 0.1 1],...
     'Width',10,'length',20,...
     'FaceColor',b0col,'EdgeColor',b0col);
 db0b=arrow([0	-0.1 1], [0 0.1 1],...
     'Width',10,'length',20,...
     'FaceColor',b0col,'EdgeColor',b0col);
 
end
%linkaxes(q1,q2)
%colorbar
hold off
axis equal


% adding labels e, i, phi, etc
fsz=18; % fontsize

ephitxt=legend('\phi','\delta E')%,'e^-','I','B_0','B_1')
set(ephitxt,'color','none','box','off',...
    'position',[0.50  0.6601    0.0930    0.0774],...
    'fontsize',fsz)
etxt1=text(1.2,0,0.5,'e^-','color',[0 0 0],'fontsize',fsz);
etxt2=text(5.2,0,0.5,'e^-','color',[0 0 0],'fontsize',fsz);
b0txt1=text(0.1,0,1.1,'B','color',[0 0 0],'fontsize',fsz);
b0txt2=text(4.1,0,1.1,'B ','color',[0 0 0],'fontsize',fsz);
b1txt1=text(0.1,0,-1.6,'\delta B','color',[0 0 0],'fontsize',fsz);
b1txt2=text(4.1,0,-1.6,'\delta B ','color',[0 0 0],'fontsize',fsz);
itxt1=text(1.2,0,-1.5,'j_{ExB}','color',[0 0 0],'fontsize',fsz);
itxt2=text(5.2,0,-1.5,'j_{ExB}','color',[0 0 0],'fontsize',fsz);

set(gca,'ylim',[-3 3],'xlim',[-2 6],'zlim',[-2 1],...
    'ycolor','white','xcolor',[1 1 1],'zcolor',[1 1 1]);
set(gcf,'color',[1 1 1])
set(gca,'CameraPosition',[2.8044  -46.0842   24.0071])
eval(['print -depsc2 db_sketchie.eps']);


%% 2d
fig=figure;
set(gcf,'PaperPositionMode','auto');
[X,Y] = meshgrid(-2:.4:2);
Z = 1.*exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,.2,.2);

[X2,Y2] = meshgrid(-2:.1:2);
Z2 = 1.*exp(-X2.^2 - Y2.^2);

r=0.5;
ang=-pi/4:0.001:pi/4;
hold on;
re=0.05;

iang=0:0.001:2*pi;
nang=size(iang,2);
ix=2+cos(iang);
iy=sin(iang);

plot3(ix,iy,ang*0,'linewidth',1.1,'color',[1 0 1])
iarrow1=arrow([vex(end) vey(end)], [vex(end) vey(end)+0.1],...
    'Length',100)
iarrow2=arrow([vex(end) vey(end)], [vex(end) vey(end)-0.1],...
    'Length',100)

%veq=quiver(sin(ang(end)),cos(ang(end)),vex*10,vey*10,1,'linewidth',1.1)
%annotation(fig,'arrow',[sin(ang(end))/4 cos(ang(end))/4],[vex vey]);
% annotation(fig,'arrow',[0.616945107398568 0.63126491646778],...
%     [0.665666666666667 0.642276422764228],'HeadLength',15,'HeadWidth',20,...
%     'HeadStyle','vback1',...
%     'LineWidth',1);

%vearrow=arrow([cos(ang(end-1)) sin(ang(end-1))], [cos(ang(end)) sin(ang(end))],...
%    'Length',100)
vearrow=arrow([vex(end) vey(end)], [vex(end)-0.1 vey(end)+0.1],...
    'Length',100)
%barrow=arrow([0 0 0.2], [0 0 0.22],...
%    'Length',10)

cont=contour(X2,Y2,Z2);
hold on
% q1=quiver(X,Y,DX,DY,0.5); set(q1,'linewidth',1)
b=quiver3(0,0,0,0,0,-0.2,1,'linewidth',1.1)
colormap hsv
hold off
axis equal
axis square
set(gca,'ylim',[-2 2],'xlim',[-2 2],'zlim',[0 0.5],...
    'ycolor','white','xcolor',[1 1 1],'zcolor',[1 1 1]);
set(gcf,'color',[1 1 1])

eval(['print -depsc2 db_sketchie.eps']);


%% 3d
figure
[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
Z = X.* exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,.1,.1);
quiver3(X,Y,Z,-DX,-DY,-DX*0,2);
hold on
surf(X,Y,Z);
colormap hsv
view(-35,45)
axis ([-2 2 -2 2 -.6 .6])
hold off
