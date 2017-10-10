for kk=1:3; h(kk)=subplot(3,1,kk); end
x0 = 0.13;
y0 = 0.13;
dy = 0.08;
yy = 0.15;

set(gcf,'position',[910    56   475   883])
set(h(1),'position',[x0 y0+2*yy+dy  0.66 0.30])
set(h(2),'position',[x0 y0+yy+dy  0.66 0.15])
set(h(3),'position',[x0 y0  0.66 0.15])

nx=500;
ny=500;
x=linspace(-4,4,nx);
y=linspace(-4,4,ny);
[X,Y] = cn.meshgrid(x,y);
Z = peaks(X,Y)';
%Z = cos(2*(X-1.5))*0.001*(1.5-0.8*abs(Y));
[Ey,Ex] = gradient(Z,x,y);

xind=fix(nx/2);xind=1:nx; 
yind=1:ny;yind=fix(nx/2)-30; 
zind=1:ny;

angles = 0:30:90;
nAngles = numel(angles);
intEk = [];
Ek = [];
ks = [tocolumn(cosd(angles)) tocolumn(sind(angles))];
for kk = 1:nAngles;       
    E = cosd(angles(kk))*Ex + sind(angles(kk))*Ey;
    intE = irf_integrate([tocolumn(1:nx) tocolumn(E(xind,yind))]);
    intEk = [intEk intE(:,2)];
    Ek = [Ek tocolumn(E(xind,yind))];
    leg{kk} = [num2str(angles(kk),'%.0f') '^o'];
end
%

set(gcf,'defaultAxesFontSize',16);
set(gcf,'defaultTextFontSize',16);
isub = 1;
if 1    
    hca=h(isub); isub=isub+1;    
    pcolor(hca,X,Y,Z); hold(hca,'on');
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'y')    
    grid(hca,'off')       
    
    plot(hca,X(xind,yind),Y(xind,yind),'k')
    quiver(hca,-2*ones(4,1),y(yind)*ones(4,1),ks(:,1),ks(:,2),1,'color',[0 0 0],'linewidth',2)
    %arrow([X(10,yind),Y(10,yind)],[1,1])
    view(hca,[0 0 1])
    irf_colormap('poynting')
    hc=colorbar('peer',hca);
    ylabel(hc,'Electrostatic potential ')
    set(hc,'ytick',1*[-1 1],'yticklabel',['-';'0'; '+'])
    axis(hca,'equal')
    set(hca,...%'xtick',[],'xticklabel',[],...
            'ztick',[],'zticklabel',[],'xlim',x([1 end]),'ylim',2*[-1 1],...
            'clim',7*[-1 1])
        
        irf_legend(hca,'spacecraft trajectory  ',[0.98 0.05])
        %annotation('arrow',[0.6 0.6],[0.62 0.73])
             annotation('arrow',[0.6 0.6],[0.62 0.65])
end
if 1
    hca=h(isub); isub=isub+1;
    plot(hca,x,Ek)
    set(hca,'ytick',[0],'ylim',max(get(hca,'ylim'))*0.7*[-1 1])
    xlabel(hca,'x')
    ylabel(hca,'Electric field')
    irf_legend(hca,leg,[1 01])
end
if 1
    hca=h(isub); isub=isub+1;
    plot(hca,x,intEk,'-',x(1:15:end),intEk(1:15:end,1),'o')
    set(hca,'ytick',[0],'ylim',max(get(hca,'ylim'))*1.1*[-1 1])
    xlabel(hca,'x')
    ylabel(hca,'Electrostatic potential')
    legs=leg;
    legs{nAngles+1} = '''B''';
    irf_legend(hca,legs,[1 1])
end 

% Labeling
abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
%legloc={[0.02,0.76],[0.02,0.8],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};

for k=1:3%nPanels
    irf_legend(h(k),[abcde(k)],[0.05 0.94],'color',colors{1})
    %grid(h(k),'off')
   % irf_legend(h(k),{'C3','C4'},[0.98 0.95]);
end
for jj = 1:3; hold(h(jj),'off'); end
irf_plot_axis_align