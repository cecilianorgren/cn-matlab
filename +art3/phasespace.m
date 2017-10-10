% plot constant energy lines of typical electron hole
close
units = irf_units;

% electrostatic and kinetic energy
u_es = @(phi)(-units.e*phi);
u_ek = @(v)(0.5*units.me*(v).^2*1e6);

% electron hole parameters
phi = @(lz,z,z0,phi0) phi0*exp(-(z-z0).^2/2/lz^2);
lz = 7;

% the box
L = 50; % km 
nz = 300;
z = linspace(0,L,nz);
vlim = 20000; nv = 200;
%vp = 1000; % km/s
vp=0;
v = linspace(-vlim-vp,vlim+vp,nv)-0*vp;
[Z,V] = meshgrid(z,v);
  
nEH = 2;
z0 = [10 30]; % km
lz = [2 3]; % km
phi0 = [400 300]; % V
vp = [11000 -2000]; 

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
nPanels = 3;
for kk=1:nPanels; h(kk)=subplot(nPanels,1,kk); end
set(gcf,'position',[95    32   472   904]);
isub=1;
if 1
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
if 1
    hca = h(isub); isub=isub+1;
    hold(hca,'on');    
    contour(hca,Z,V,U_ek);%,-3:0.5:20); 
   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x [km]')
    ylabel(hca,'v_x [km/s]')
    title(hca,'U_{ek}')
    colorbar('peer',hca)
    hold(hca,'off');    
end
if 1
    hca = h(isub); isub=isub+1;
    hold(hca,'on');    
    contour(hca,Z,V,U);%,-3:0.5:20); 
   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x [km]')
    ylabel(hca,'v_x [km/s]')
    title(hca,'U_{ek}+U_{es}')
    colorbar('peer',hca)
    hold(hca,'off');    
end
%%
    title(hca,'U_{ek}+U_{es}')
    %colorbar('peer',hca)
    set(hca,'xtick',L*[0:nL],'xticklabel',{0 1 2 3 4}); xlabel(hca,'x/\lambda')
    yticks = get(hca,'ytick');
    yticks = min(yticks):2:max(yticks);
    yticklabels = num2str(yticks');
    set(hca,'ytick',yticks,'yticklabel',yticklabels); ylabel(hca,'v_x/v_p');
    set(hca,'ylim',max(abs(get(hca,'ylim')))*[-1 1])
    %set(gcf,'position',[930   563   601   250])
    box(hca,'on')
    
end
if nPanels == 2
    hca = h(isub); isub=isub+1;
    plot(x,phi(k,x,o,t))
    %contour(hca,X,V,U,-1:0.5:5); 
    %pcolor(hca,X,V,U); 
    %view(hca,[0 0 1]); 
    %shading(hca,'flat');
    %xlabel(hca,'x')
    %ylabel(hca,'v_x')
    title(hca,'')
    
    %set(hca,'
end
%h(2)=subplot(3,1,2); surf(X,V,U_es); view([0 0 1]); shading flat;

%h(3)=subplot(3,1,3); surf(X,V,U_ek); view([0 0 1]); shading flat;


