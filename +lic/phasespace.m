% plot constant energy lines of typical electron hole
close
units = irf_units;

% electrostatic and kinetic energy
u_es = @(phi)(-units.e*phi);
u_ek = @(v)(0.5*units.me*(v).^2*1e6);

% electron hole parameters
phi = @(lz,z,z0,phi0) phi0*exp(-(z-z0).^2/2/lz^2);
lz = 5;

% the box
L = 25; % km 
nz = 300;
z = linspace(0,L,nz);
vlim = 30000; nv = 200;
%vp = 1000; % km/s
vp=0;
v = linspace(-vlim-vp,vlim+vp,nv)-0*vp;
[Z,V] = meshgrid(z,v);
  
nEH = 1;
z0 = [10 30]; % km
lz = [2 3]; % km
phi0 = [500 300]; % V
vp = [500 -2000]; 

U_es=0;
U_ek=0;
U = 0;
for ii=1:nEH
    vind = find(V>vp(ii),1,'first')
    U_es = U_es + u_es(phi(lz(ii),Z,z0(ii),phi0(ii)));
    U = U + u_es(phi(lz(ii),Z,z0(ii),phi0(ii))) + u_ek(V-vp(ii));    
    %U_ek = 
end
U_ek = u_ek(V);     
%U = (+ U_es + U_ek);
Udiff = (- U_es + U_ek);

% set up figure
nPanels = 1;
for kk=1:nPanels; h(kk)=subplot(nPanels,1,kk); end
%set(gcf,'position',[95    32   472   904]);
isub=1;
if 0
    hca = h(isub); isub=isub+1;
    hold(hca,'on');    
    contour(hca,Z,V,U_es);   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x [km]')
    ylabel(hca,'v_x [km/s]')
    title(hca,'U_{es}')
    colorbar('peer',hca)
    hold(hca,'off');    
end
if 0
    hca = h(isub); isub=isub+1;
    hold(hca,'on');    
    contour(hca,Z,V,U_ek);%,-3:0.5:20); 
   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x [km]')
    ylabel(hca,'v_{x} [km/s]')
    title(hca,'U_{ek}')
    colorbar('peer',hca)
    hold(hca,'off');    
end
if 1
    hca = h(isub); isub=isub+1;
    hold(hca,'on');    
    contour(hca,Z,V*1e-3,U*1e16,30);%,-3:0.5:20); 
   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x [km]')
    ylabel(hca,'v_x [10^3 km/s]')
    title(hca,'U_{ek}+U_{es}')
    %ch = colorbar('peer',hca);
    %set(ch,'ylabel','J')
    hold(hca,'off');  
    box(hca,'on')
end


