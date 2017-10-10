% fig 3 in article directly combined



% load electric field
cd /Users/Cecilia/Data/BM/20070831
load fac_dxdydz
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields

if ~exist('ind','var'); ind = 5; end
ind = 5; % Has nice mono.
tint = tint{highQuality(ind)};

if 0 % the electric field
epar4 = irf_tlim(EparAC4,tint);
epar3 = irf_tlim(EparAC3,tint);
eper4 = irf_tlim(Epar4,tint);
eper3 = irf_tlim(Epar3,tint);
end
if 1 % the coordinate system
dx=cn.mean(irf_tlim(dp1,tint),1); % SP
dy=cn.mean(irf_tlim(dp2,tint),1); % BxSP
dz=cn.mean(irf_tlim(dp3,tint),1); % B
end

% model from Chen2004: Bernstein-Greene-Kruskal solitary waves in...
phi = @(phi0,r,z,lz,lr) phi0*exp(-z.^2/2/lz.^2-r.^2/2/lr.^2);
% deduced fit parameters
lz = 5; % km
lr = 12; % km
phi0 = 500; % V
%%
% Plot spacecraft spin plane and E-components
% Sketch, fig 3a in paper
subplot(1,2,1)

spangle = 30;
rspangle = 0.2;
spanglee = linspace(0,spangle,30);

angle = linspace(0,2*pi,100);
r = 1;
Ereal = 1;
sca = 1.05;
arr = 1.1;

fsize = 22;
% spin plane
plot3(r*cos(angle),r*sin(angle),r*angle*0,'k'); hold on;
text(-r*0.2,r*0.7,r*0,{'spin','plane'},'fontsize',fsize)

% b-field
bx = cosd(spangle)*[-0.7 1]*1.2;
by = [0 0];
bz = sind(spangle)*[-0.7 1]*1.2;
plot3(r*bx,r*by,r*bz);
text(-r*0.3,r*0.2,-r*0.1,'B','fontsize',fsize)

% angle
plot3(r*rspangle*cosd(spanglee),spanglee*0,r*rspangle*sind(spanglee),'color',[0 0 0])
text(r*0.3,r*0,r*0.1,'\theta','fontsize',fsize+4)

% e-par
quiver3(0,0,0,r*Ereal*cosd(spangle),r*0,r*Ereal*sind(spangle),sca,'r','linewidth',2)
text(r*(Ereal*cosd(spangle)-0.4),0,r*Ereal*sind(spangle)-0.0,'E_{||}','fontsize',fsize)
%arrow3([0 0 0],arr*[r*Ereal*cosd(spangle) 0 r*Ereal*sind(spangle)]+[0,0,0])

% e-par,sp
quiver(0,0,r*Ereal*cosd(spangle),r*0,sca,'g','linewidth',2)
text(r*Ereal*cosd(spangle)*0.65,-r*0.2,'E_{||,meas}','fontsize',fsize)
%arrow3([0 0 0],arr*[r*Ereal*cosd(spangle) 0 0]+[0,0,0])

% e-per,sp
quiver(0,0,0,-r*Ereal*cosd(spangle),sca,'g','linewidth',2)
%arrow3([0 0 0],arr*[0 -r*Ereal*cosd(spangle) 0]+[0,0,0])
%arrow3([x4 y4 0],[x4 y4 ]+[5,0,0],'length',10,'baseangle',70)
text(r*0.1,-r*0.7,'E_{\perp,meas}','fontsize',fsize)

% projection line
line(arr*r*[Ereal*cosd(spangle) Ereal*cosd(spangle)],[0 0],arr*r*[0 Ereal*sind(spangle)],'linestyle','--','linewidth',2,'color',[0 0 0])


if 1
    text(-r*1.25,r*0,r*0.97,'a)','fontsize',28)
end

set(gca,'xlim',[-1 1],'ylim',[-1 1])
zoom(gca,1.35)



axis equal
axis square
axis off

hold off;



% Figure 3b
subplot(1,2,2)

x=linspace(-lr*1.55,2*lr,200);
y=linspace(-3.8*lr,2*lr,200);
[X,Y]=meshgrid(x,y);

fs=20;

R = sqrt(X.^2+Y.^2);
levels = phi0*[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 ];
[C,h] = contour(X,Y,phi(phi0,R,0,lz,lr),levels);
cl = clabel(C,h,levels(1:2:end),'labelspacing',1000,'fontsize',fs);
for kk=1:numel(cl);
    str = get(cl(kk),'String');
    str = [str ' V'];
    set(cl(kk),'String',str);
    set(cl(kk),'BackgroundColor','white');
end
xlabel('in spin plane (\perp B)')
ylabel('out of spin plane (\perp B)')
xlabel('[km]','fontsize',fs)
ylabel('[km]','fontsize',fs)
hold on;
% add sc position
% dp = p3-p4

x3 = 10; x4 = x3-dx; % in spin plane
y3 = 10; y4 = y3-dy; % out of spin plane
plot(x3,y3,'og',x4,y4,'vb','markersize',20,'linewidth',2)
text(x3-5,y3,'C3','fontsize',fs+2)
text(x4-5,y4,'C4','fontsize',fs+2)
axis equal


arrow([x3 y3],[x3 y3]+[5,0],'length',10,'baseangle',70)
arrow([x4 y4],[x4 y4]+[5,0],'length',10,'baseangle',70)
text(x4+4,y4-2,'E_{\perp,meas}','fontsize',fs+2)
text(x3+4,y3-2,'E_{\perp,meas}','fontsize',fs+2)


if 1 % add b) label
    text(-16,13,'b)','fontsize',28)
end
    
%axis square
set(gca,'ylim',[-20 15],'xlim',[-9 20])
set(gca,'fontsize',fs)
hold off;



if 0
% make electric field to compare
zlim = 30; % km
nz = 100;
z = linspace(-zlim,zlim,nz);
theta = linspace(0,359,360);
rlim = 30; nr = 50;
r = linspace(0,rlim,nr);
x = r'*cosd(theta);
y = r'*sind(theta);

[X Y Z] = meshgrid(x,y,z);
% define trajectory, radius at which eh pass
rt = 5; % km
end