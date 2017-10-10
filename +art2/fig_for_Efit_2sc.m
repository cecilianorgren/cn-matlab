
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


%% first version of plot, seen from B direction with equipotential lines
% make plot

set(0,'defaultTextFontSize',18);
set(0,'defaultAxesFontSize',18);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');
        
x=linspace(-lr*1.60,2*lr,200);
y=linspace(-2*lr,2*lr,200);
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
    text(-15,13,'b)','fontsize',fs+1)
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

%% second version of plot, box seen from angle with equipotential lines, but only in the directions we can measure
% make plot
x=linspace(-lr*0.60,2*lr,20);
y=linspace(-2*lr,2*lr,20);
z=linspace(-3*lz,3*lz,20);
[Xz,Yz]=meshgrid(x,y);
[Yx,Zx]=meshgrid(y,z);
[Xy,Zy]=meshgrid(x,z);
[X,Y,Z]=meshgrid(x,y,z);
R = sqrt(X.^2+Y.^2);


% add sc position
% dp = p3-p4
x3 = 10; x4 = x3-dx; % in spin plane
y3 = 10; y4 = y3-dy; % out of spin plane
z3 = 10; z4 = z3-dz; % along b
plot3(x3,y3,z3,'og',x4,y4,z4,'vb','markersize',20,'linewidth',2)
hold on;

ax=gca;
a = get(ax, 'XLim'); xpos = a(1);
a = get(ax, 'YLim'); ypos = a(1);
a = get(ax, 'ZLim'); zpos = a(1);

xlabel('x')
ylabel('y')
zlabel('z')
%
if 1 % projection on bottom z plane
    Rz = sqrt(Xz.^2+Yz.^2);    
    PHIz = phi(phi0,Rz,0,lz,lr);
    [~,hh]=contour3(ax,Xz,Yz,PHIz,levels);
    % size zpos to match the data
    for i = 1 : length(hh)
        zz = get(hh(i), 'ZData');
        set(hh(i), 'ZData', zpos * ones(size(zz)));
    end
end
if 1 % projection on bottom x plane    
    Rx = sqrt(0*Xx.^2+Yx.^2);    
    PHIx = phi(phi0,Rx,Zx,lz,lr);
    [~,hh]=contour3(ax,PHIx,Yx,Zx,levels);

    % size zpos to match the data
    for i = 1 : length(hh)
        xx = get(hh(i), 'xData');
        set(hh(i), 'xData', xpos * ones(size(xx)));
    end
end
if 1 % projection on bottom y plane    
    Ry = sqrt(Xy.^2+0*Yz.^2);
    PHIy = phi(phi0,Xy,Zy,lz,lr);
    [~,hh]=contour3(ax,Xy,PHIy,Zy);

    % size zpos to match the data
    for i = 1 : length(hh)
        yy = get(hh(i), 'yData');
        set(hh(i), 'yData', ypos * ones(size(yy)));
    end
end
hold(ax,'off')

%%
Rz = sqrt(Xz.^2+Yz.^2);
Ry = sqrt(Xy.^2+Yy.^2);
R = sqrt(X.^2+Y.^2);
levels = phi0*[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 ];
contour3(X,Y,phi(phi0,R,0,lz,lr),levels);
PHIz = phi(phi0,Rz,0,lz,lr);
PHIy = phi(phi0,Ry,0,lz,lr);
%hz = surfc(Xz,Yz,PHIz);
%hy = surfc(Xy,Yy,PHIy);
clabel(C,h,levels(1:2:end),'labelspacing',1000);
xlabel('in sp (\perp B)')
ylabel('out of sp (\perp B)')

delete(findobj(gcf,'type','surf'))
%%
% Plot surface
        hs = surf(cax, args{:});
        
        hold(cax, 'on');
        
        a = get(cax, 'ZLim');
        
        % Always put contour below the plot.
        zpos = a(1);
        
        % Get D contour data
        [~, hh] = contour3(cax, x, y, z);
        
        % size zpos to match the data
        for i = 1 : length(hh)
            zz = get(hh(i), 'ZData');
            set(hh(i), 'ZData', zpos * ones(size(zz)));
        end
%%
axis equal
%axis square
set(gca,'ylim',[-20 20],'xlim',[-5 20])
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