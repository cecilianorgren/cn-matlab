kmin = 0;
kmax = 1; 
nk=100;
kper = linspace(kmin,kmax,nk);
kpar = linspace(kmin,kmax,nk);

[KPER,KPAR] = meshgrid(kper,kpar);
N = 1;
om = @(kper,kpar) kper*N./sqrt(kper.^2+kpar.^2);

OMEGA = om(KPER,KPAR);

contour(KPER,KPAR,OMEGA',20)
title({'Lines of constant omega (given \theta)',...
       'Lines coincides with phase velocity lines || k'})
xlabel('k_{\perp}')
ylabel('k_{||}')

%
hold on;
ind=2:9:nk;
[vgpar,vgper] = gradient(OMEGA(ind,ind));
quiver(KPER(ind,ind),KPAR(ind,ind),KPER(ind,ind),KPAR(ind,ind))
quiver(KPER(ind,ind),KPAR(ind,ind),vgper',vgpar')

hold off
legend(' - omega(k_{\perp},k_{||})','v_p','v_g')

%% surface plot
kmin = 0;
kmax = 1; 
nk=50;
kper = linspace(kmin,kmax,nk);
kpar = linspace(kmin,kmax,nk);

[KPER,KPAR] = meshgrid(kper,kpar);
N = 1;
om = @(kper,kpar) kper*N./sqrt(kper.^2+kpar.^2);

OMEGA = om(KPER,KPAR);

surf(KPER,KPAR,OMEGA')
title({'Dispersion surface, \omega = \omega(k_{\perp},k_{||})',...
       'The direction of v_p=\omega/k is given by a line from the origin out to a given (k_{\perp},k_{||})',...
       'The direction of v_g=\partial\omega/\partial k is given by the slope of the surface.'})
xlabel('k_{\perp}')
ylabel('k_{||}')
zlabel('\omega')
