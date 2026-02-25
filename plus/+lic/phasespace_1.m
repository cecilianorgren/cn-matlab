%

units = irf_units;
u_es = @(phi,e)(e*(phi));
u_ek = @(m,v,vp)(0.5*m*(v-vp).^2);

L = 4; nx = 300;
nL = 3;
k = 2*pi/L;
o = 0.5*pi;
vp = o/k;
x = linspace(0,nL*L,nx);
l = x/L;

vlim = 3; nv = 200;
v = linspace(-vlim-vp,vlim+vp,nv)-0*vp;

%T = 1; o = 2*2*pi/T; nt = 30;
%t = linspace(0,0.5*T,nt);

[X,V] = meshgrid(x,v);
Xl= X/L;
%phi = @(k,x,o,t) cos(k*x-o*t+t*o+pi)+0*1;%.*(1+0.1*x);
phi = @(k,x,o,t) cos(k*x)+0*1;%.*(1+0.1*x);
e = -1;
m = 1;
t = 0;

V=V;

nPanels = 1;
for kk=1:nPanels; h(kk)=subplot(nPanels,1,kk); end
    
U_es = u_es(phi(k,X,o,t),e) + 0*u_es(phi(-k,X,o,t),e);
U_ek = u_ek(m,V,vp) + 0*u_ek(m,V,vp); 
U = (+ U_es + U_ek);
Udiff = (- U_es + U_ek);

isub=1;
if 1
    hca = h(isub); isub=isub+1;
    hold(hca,'on');
    %pcolor(hca,X,V,Udiff); 
    contour(hca,X,V,U,-3:0.5:20); 
   
    view(hca,[0 0 1]); 
    shading(hca,'flat');
    xlabel(hca,'x')
    ylabel(hca,'v_x')
    %title(hca,'U_{ek}+U_{es}')
    title(hca,'\epsilon + q\phi','interpreter','tex')
    %colorbar('peer',hca)
    set(hca,'xtick',L*[0:nL],'xticklabel',{0 1 2 3 4}); xlabel(hca,'x/\lambda')
    yticks = get(hca,'ytick');
    yticks = min(yticks):2:max(yticks);
    yticklabels = num2str(yticks');
    set(hca,'ytick',yticks,'yticklabel',yticklabels); ylabel(hca,'v_x/v_{ph}');
    set(hca,'ylim',max(abs(get(hca,'ylim')))*[-1 1])
    %set(gcf,'position',[930   563   601   250])
    box(hca,'on')
    hca.FontSize = 14;
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


