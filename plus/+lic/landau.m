f = @(n,v,vt) (n/(2*vt)^(3/2)).*exp(-v.^2/(vt.^2));
vt = cn_eV2v(1600,'eV');

nv = 1000;
v = linspace(-3,3,nv)*vt;
accf = zeros(nv,1);
n = 0.04*1e6; 

vd=2*vt;
vtb=0.4*vt;
nb=0.1*n;

set(gcf,'defaultAxesFontSize',16);
subplot(1,1,1)
ftot=f(n,v,vt)+f(nb,v-vd,vtb);
plot(v/vt,ftot,'k','linewidth',1); hold on;
title('f(v)','fontsize',18)
plot(v/vt,v*0,'k','linewidth',1)
text(0,-0.08*ftot(500),'v','fontsize',20)
set(gca,'xtick',[],'ytick',[])
%xlabel('v','fontsize',18)
%
vsta=find(v>1*vt,1,'first');
vsto=find(v<1.1*vt,1,'last');
vind=vsta:vsto;
nvind=numel(vind);
X = [flipdim(v(vind),2) v(vind) ];
Y = [flipdim(ftot(vind),2) zeros(1,nvind) ];
patch(X/vt,Y,[1 0.7 0.4])
text(v(vsto)/vt*0.9,ftot(vsto)*1.6,'v_{ph,1}','fontsize',18)
%vmid=vind(fix(end/2))+0;
%plot([v(vmid) v(vmid)]/vt,[0 ftot(vmid)],'k')

vsta=find(v>1.1*vt,1,'first');
vsto=find(v<1.2*vt,1,'last');
vind=vsta:vsto;
nvind=numel(vind);
X = [flipdim(v(vind),2) v(vind) ];
Y = [flipdim(ftot(vind),2) zeros(1,nvind) ];
patch(X/vt,Y,[1 0.5 0.7])
%vmid=vind(fix(end/2))+0;
%plot([v(vmid) v(vmid)]/vt,[0 ftot(vmid)],'k')


vsta=find(v>1.65*vt,1,'first');
vsto=find(v<1.75*vt,1,'last');
vind=vsta:vsto;
nvind=numel(vind);
X = [flipdim(v(vind),2) v(vind) ];
Y = [flipdim(ftot(vind),2) zeros(1,nvind) ];
patch(X/vt,Y,[1 0.7 0.4])
vmid=vind(fix(end/2))+0;
%plot([v(vmid) v(vmid)]/vt,[0 ftot(vmid)],'k')
text(v(vsto)/vt*0.9,ftot(vsto)*1.7,'v_{ph,2}','fontsize',18)
vsta=find(v>1.75*vt,1,'first');
vsto=find(v<1.85*vt,1,'last');
vind=vsta:vsto;
nvind=numel(vind);
X = [flipdim(v(vind),2) v(vind) ];
Y = [flipdim(ftot(vind),2) zeros(1,nvind) ];
patch(X/vt,Y,[1 0.5 0.7])
vmid=vind(fix(end/2))+0;

if 0
plot(v/vt,ftot,'k','linewidth',2); hold on;
title('f(v)','fontsize',18)
plot(v/vt,v*0,'k','linewidth',2)
text(0,-0.08*ftot(500),'v','fontsize',20)
end
hold off


%%
subplot(2,1,2)
plot(v/vt,accf)
title('accumulated f(v)')

%% 2d

f = @(n,vx,vy,vt) (n/(2*vt)^(3/2)).*exp(-(vx.^2+vy.^2)/(vt.^2));
vt = cn_eV2v(1600,'eV');

nv = 30;
vx = linspace(0,2,nv)*vt;
vy = linspace(0,2,nv)*vt;
accf = zeros(nv,nv);
n = 0.04*1e6; 

vd=0.4*vt;
vtb=0.2*vt;
nb=0.00*n;
for kk = 1:nv
    for pp = 1:nv
        lx=sort([0,vx(kk)]);
        ly=sort([0,vy(pp)]);
    accf(kk,pp) = integral2(@(vx,vy)f(n,vx,vy,vt),lx(1),lx(2),ly(1),ly(2));
    end
end
%%
[VX,VY]=meshgrid(vx,vy);
subplot(2,1,1)
pcolor(VX,VY,accf')
shading flat
colorbar
axis square

subplot(2,1,2)
pcolor(VX,VY,f(n,VX,VY,vt)')
title('f(v)')
axis square
%%



plot(v/vt,f(n,v,vt)+f(nb,v-vd,vtb))

subplot(2,1,2)
plot(v/vt,accf)
title('accumulated f(v)')

%
