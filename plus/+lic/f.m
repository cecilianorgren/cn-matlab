f = @(n,v,vt) (n/(2*vt)^(3/2)).*exp(-v.^2/(vt.^2));
vt = cn_eV2v(1600,'eV');

nv = 300;
v = linspace(-3,3,nv)*vt;
accf = zeros(nv,1);
n = 0.04*1e6; 

vd=0.4*vt;
vtb=0.2*vt;
nb=0.00*n;
for kk = 1:nv
    accf(kk) = integral(@(v) f(n,v,vt)+f(nb,v-vd,vtb),v(1),v(kk));
end

subplot(2,1,1)
plot(v/vt,f(n,v,vt)+f(nb,v-vd,vtb))
title('f(v)')
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

%% For licentiate presentation
f = @(n,v,vd,vt) (n^(1)/((pi^(1/2)*vt)^(1))).*exp(-(v-vd).^2/(1*vt.^2));

units = irf_units;

if 0
    vte1 = cn_eV2v(1600,'eV');
    vte2 = cn_eV2v(0.95,'eV');
    vti = cn_eV2v(1600,'eV')*sqrt(units.me/units.mp); 

    nv = 500;
    v = linspace(-3,3,nv)*vte1;
    n = 0.04*1e6; 
    R = 1;
    n2 = n*R;
    n1 = n*(1-R);

    vd1=0*vt;
    vd2=0.4*vte1;
    vdi=0;
    fiscale = 1;
else
    vte1 = cn_eV2v(1600,'eV');
    vte2 = cn_eV2v(60,'eV');
    vti = cn_eV2v(1600,'eV')*sqrt(units.me/units.mp); 

    nv = 500;
    v = linspace(-5,5,nv)*vte1;
    n = 0.04*1e6; 
    R = 0.0;
    n2 = n*R;
    n1 = n*(1-R);

    vd1=1.9*vte1;
    vd2=.4*vte1;
    vdi=0;
    fiscale = 0.05;
end
if 1
    h = plot(v/vte1,f(n1,v,vd1,vte1)+f(n2,v,vd2,vte2),...
        v/vte1,fiscale*f(n,v,vdi,vti),...
        v/vte1,f(n2,v,vd2,vte2),...
        v/vte1,f(n1,v,vd1,vte1),...
        'linewidth',2);
    %t('f(v_{||})')
    xlabel('v_{||}/v_{te}')
    irf_legend({' ','f_i'},[0.38 0.8])
    irf_legend({'f_e',' '},[0.7 0.5])
    set(gca,'ytick',[])
    set(gcf,'position',[103   485   423   225]);
else
    for kk = 1:nv    
       accf(kk,1) = integral(@(v) f(n1,v,vd1,vte1),v(1),v(kk));
       accf(kk,2) = integral(@(v) f(n2,v,vd2,vte2),v(1),v(kk));
    end
    plot(v/vte1,accf)
end