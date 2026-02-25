%
close
units = irf_units;
u_es = @(phi,e)(e*phi);
u_ek = @(m,v,vp)(0.5*m*(v-vp).^2);

vlim = 3; nv = 200;
v = linspace(-vlim,vlim,nv);

L = 4; nx = 300;
k = 2*pi/L;
x = linspace(0,3*L,nx);

T = 1; o = 2*2*pi/T; nt = 30;
t = linspace(0,0.5*T,nt);

[X,V] = meshgrid(x,v);

phi = @(k,x,o,t) cos(k*x-o*t+t*o);

e = 1;
m = 1;
vp = o/k;
V=V+0*vp;

for ii = 1%:nt
    U_es = u_es(phi(k,X,o,t(ii)),e) + ...
           0*u_es(phi(-k,X,o,t(ii)),e);
       
    U_ek = u_ek(m,V,-vp*0) +...
           0*u_ek(m,V,vp); 
    
       U = U_es + U_ek;
    if 0
        contour(X,V,U,12,'linewidth',1.1);  shading flat; %colorbar
        hold on; [c,b] = contour(X,V,U,[1,1],'linewidth',1.1); %clabel(c,b);
        axis equal; set(gca,'ylim',[-vlim vlim],'xlim',[x([1 end])])
        xlabel(gca,'x')
        ylabel(gca,'v_x')
    else
        h(1)=subplot(3,1,1); surf(X,V,U); view([0 0 1]); shading flat;
        h(2)=subplot(3,1,2); surf(X,V,U_es); view([0 0 1]); shading flat;
        cb2 = colorbar; caxis(2*[-1 1])
        h(3)=subplot(3,1,3); surf(X,V,U_ek); view([0 0 1]); shading flat;
        pause(0.1)
        for jj =1:3
            xlabel(h(jj),'x')
            ylabel(h(jj),'v_x')
        end
    end
    
end
