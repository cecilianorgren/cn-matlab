%% make model electric field
e = 1.6022e-19;
me = 9.1094e-31;
n = 0.08*1e6; % m^-3
mu0 = 1.257e-6;
B0 = 20*1e-9; % T

% cosine wave
phi0 = 300; % V
phi = @(k,x,w,t) phi0*cos(k*x-w*t);
E = @(k,x,w,t,phase) -k*phi0*sin(k*x-w*t+phase);
v_par = @(k,x,w,t,phase) e*E(k,x,w,t,phase)/w/me; % m/s
% j_par = -n*e*v_par;

% set up system
vph = -1400*1e3; % m/s
rhoe = 10e3; % m
lam = 2*pi*rhoe; % m
k = 2*pi/lam; % m^-1
w = vph*k;
T = 2*pi/w;
t = linspace(0,3*T,150);
x = 0;

prop_angle = 89*pi/180;
kper = k*sin(prop_angle);
kpar = k*cos(prop_angle);

% calculate current from velocity fields
% tool.jXB
phi2 = @(k,x,y,w,t) phi0*cos(k*x-w*t).*sin(k*y);
x = t*vph;
y = linspace(0,lam,100);
z = linspace(-mean(x),+mean(x),20);
[sX,sY,sZ] = meshgrid(x,y,z);
[fX,fY] = meshgrid(x,y);
PHI = phi2(k,0,fY,w,fX/vph);
sJZ = 0.14e-9*PHI; % A, jpar=0.14e-9*phi, scaling factor, CAN NOT USE THIS since it should depend on k_angle
sJZ = e*n*e*kpar/w/me*PHI;
[BX BY] = tool.j2B(repmat(sJZ,1,1,size(sX,3)),sX,sY,sZ,fX,fY);
BZ = PHI*e*n*mu0/B0;
% [BX BY] = tool.j2B(sJZ,sX,sY,fX,fY);
maxB = max([max(max(BX)) max(max(BY)) max(max(BZ))]);
%%
figure(2)
npl=4;
for kk=1:npl; h(kk)=subplot(npl,1,kk); end
isub = 1;

if 1
    hca=h(isub); isub=isub+1;
    hold(hca,'on')
    ply = 1:3:numel(t); % plot step
    plx = 1:3:numel(y); % plot step
    contour(hca,fX*1e-3,fY*1e-3,BZ)
    quiver(hca,fX(plx,ply)*1e-3,fY(plx,ply)*1e-3,BX(plx,ply),BY(plx,ply))
    ylabel(hca,'y [km]')
    xlabel(hca,'x (=vph*t) [km]')
    title(h(1),['LHDW B reconstruction, \theta = ' num2str(prop_angle*180/pi) '^o',...
                ', v_{ph}=' num2str(vph*1e-3) ' km/s'])
    box(hca,'on')
    axis equal
    set(hca,'ylim',y([1 end])*1e-3,'xlim',vph*t([1 end])*1e-3)
end
if 1
    hca=h(isub); isub=isub+1;
    yind = 13;
    plot(hca,x*1e-3,squeeze(BX(yind,:))*1e9,...
             x*1e-3,squeeze(BY(yind,:))*1e9,...
             x*1e-3,squeeze(BZ(yind,:))*1e9)
    plot(h(1),x*1e-3,x*0+y(yind)*1e-3,'k')     
    ylabel(hca,['B at y=' num2str(y(yind)*1e-3,'%.0f') ' [nT]'])
    irf_legend(hca,{'X','Y','Z'},[1 0.95])
    set(hca,'xlim',x([1 end])*1e-3,'ylim',maxB*1e9*[-1 1])%,'xlim',vph*t([1 end])*1e-3)
    Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) tocolumn(squeeze(BZ(yind,:))*1e9)];
    [~,l,v]=irf_minvar(Bt);
    ang2B=acosd(abs([0 0 1]*v(3,:)'));
    ang2k=acosd(abs([1 0 0]*v(3,:)'));
    title(hca,['v1=[' num2str(v(1,1),'%.2f') ' ' num2str(v(1,2),'%.2f') ' ' num2str(v(1,3),'%.2f') ']',...
        ', v2=[' num2str(v(2,1),'%.2f') ' ' num2str(v(2,2),'%.2f') ' ' num2str(v(2,3),'%.2f') ']',...
        ', v3=[' num2str(v(3,1),'%.2f') ' ' num2str(v(3,2),'%.2f') ' ' num2str(v(3,3),'%.2f') ']',...
        ', L1/L2=' num2str(l(1)/l(2),'%1.1e\n') ', L2/L3=' num2str(l(2)/l(3),'%1.1e\n'),...  
        ', \theta_B=' num2str(ang2B,'%.0f') '^o, \theta_x=' num2str(ang2k,'%.0f') '^o'])
    end
if 1
    hca=h(isub); isub=isub+1;
    yind = 27;
    plot(hca,x*1e-3,squeeze(BX(yind,:))*1e9,...
             x*1e-3,squeeze(BY(yind,:))*1e9,...
             x*1e-3,squeeze(BZ(yind,:))*1e9)
    plot(h(1),x*1e-3,x*0+y(yind)*1e-3,'k')          
    ylabel(hca,['B at y=' num2str(y(yind)*1e-3,'%.0f') ' [nT]'])
    irf_legend(hca,{'X','Y','Z'},[1 0.95])
    set(hca,'xlim',x([1 end])*1e-3,'ylim',maxB*1e9*[-1 1])%,'xlim',vph*t([1 end])*1e-3)
    Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) tocolumn(squeeze(BZ(yind,:))*1e9)];
    [~,l,v]=irf_minvar(Bt);
    ang2B=acosd(abs([0 0 1]*v(3,:)'));
    ang2k=acosd(abs([1 0 0]*v(3,:)'));
    title(hca,['v1=[' num2str(v(1,1),'%.2f') ' ' num2str(v(1,2),'%.2f') ' ' num2str(v(1,3),'%.2f') ']',...
        ', v2=[' num2str(v(2,1),'%.2f') ' ' num2str(v(2,2),'%.2f') ' ' num2str(v(2,3),'%.2f') ']',...
        ', v3=[' num2str(v(3,1),'%.2f') ' ' num2str(v(3,2),'%.2f') ' ' num2str(v(3,3),'%.2f') ']',...
        ', L1/L2=' num2str(l(1)/l(2),'%1.1e\n') ', L2/L3=' num2str(l(2)/l(3),'%1.1e\n'),...  
        ', \theta_B=' num2str(ang2B,'%.0f') '^o, \theta_x=' num2str(ang2k,'%.0f') '^o'])
end
if 1
    hca=h(isub); isub=isub+1;
    yind = 40;
    plot(hca,x*1e-3,squeeze(BX(yind,:))*1e9,...
             x*1e-3,squeeze(BY(yind,:))*1e9,...
             x*1e-3,squeeze(BZ(yind,:))*1e9)
    plot(h(1),x*1e-3,x*0+y(yind)*1e-3,'k')     
    ylabel(hca,['B at y=' num2str(y(yind)*1e-3,'%.0f') ' [nT]'])
    irf_legend(hca,{'X','Y','Z'},[1 0.95])
    set(hca,'xlim',x([1 end])*1e-3,'ylim',maxB*1e9*[-1 1])%,'xlim',vph*t([1 end])*1e-3)
    Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) 0*tocolumn(squeeze(BZ(yind,:))*1e9)];
    [~,l,v]=irf_minvar(Bt); 
    
    ang2B=acosd(abs([0 0 1]*v(3,:)'));
    ang2k=acosd(abs([1 0 0]*v(3,:)'));
    title(hca,['v1=[' num2str(v(1,1),'%.2f') ' ' num2str(v(1,2),'%.2f') ' ' num2str(v(1,3),'%.2f') ']',...
        ', v2=[' num2str(v(2,1),'%.2f') ' ' num2str(v(2,2),'%.2f') ' ' num2str(v(2,3),'%.2f') ']',...
        ', v3=[' num2str(v(3,1),'%.2f') ' ' num2str(v(3,2),'%.2f') ' ' num2str(v(3,3),'%.2f') ']',...
        ', L1/L2=' num2str(l(1)/l(2),'%1.1e\n') ', L2/L3=' num2str(l(2)/l(3),'%1.1e\n'),...  
        ', \theta_B=' num2str(ang2B,'%.0f') '^o, \theta_x=' num2str(ang2k,'%.0f') '^o'])
end

hold(h(1),'off')
%% make quiver plot of magnetic field to understand minvar change
figure(2)
npl=4;
for kk=1:npl; h(kk)=subplot(1,npl,kk); end
isub = 1;
yinds = [13 20 27 51];
for kk=1:numel(yinds);
if 1
    hca=h(isub); isub=isub+1;
    yind = yinds(kk);40;    
    pl = 1:3:numel(t); % plot step  
    nq = numel(pl);
    orig = zeros(1,nq);
    quiver3(hca,orig,orig,orig,...
                squeeze(BX(yind,pl))*1e9,...
                squeeze(BY(yind,pl))*1e9,...
                squeeze(BZ(yind,pl))*1e9)        
    xlabel(hca,'x [nT]');ylabel(hca,'y [nT]');zlabel(hca,'z [nT]')
    irf_legend(hca,{['B_{y=' num2str(y(yind)/lam,'%.2f') '} [nT]']},[1 0.95])
    %set(hca,'xlim',x([1 end])*1e-3,'ylim',maxB*1e9*[-1 1])%,'xlim',vph*t([1 end])*1e-3)
   % axis(hca,'equal')
    view(hca,[0 1 0])
    Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) tocolumn(squeeze(BZ(yind,:))*1e9)];
    [~,l,v]=irf_minvar(Bt);
    ang2B=acosd(abs([0 0 1]*v(3,:)'));
    ang2k=acosd(abs([1 0 0]*v(3,:)'));
    titlestr = {['v1=[' num2str(v(1,1),'%.2f') ' ' num2str(v(1,2),'%.2f') ' ' num2str(v(1,3),'%.2f') ']',...
                ', v2=[' num2str(v(2,1),'%.2f') ' ' num2str(v(2,2),'%.2f') ' ' num2str(v(2,3),'%.2f') ']',...
                ', v3=[' num2str(v(3,1),'%.2f') ' ' num2str(v(3,2),'%.2f') ' ' num2str(v(3,3),'%.2f') ']',...
                ', L1/L2=' num2str(l(1)/l(2),'%1.1e\n') ', L2/L3=' num2str(l(2)/l(3),'%1.1e\n'),...  
                ', \theta_B=' num2str(ang2B,'%.0f') '^o, \theta_x=' num2str(ang2k,'%.0f') '^o']};
    titlestr = {['v1=[' num2str(v(1,1),'%.2f') ' ' num2str(v(1,2),'%.2f') ' ' num2str(v(1,3),'%.2f') ']'],...
                ['v2=[' num2str(v(2,1),'%.2f') ' ' num2str(v(2,2),'%.2f') ' ' num2str(v(2,3),'%.2f') ']'],...
                ['v3=[' num2str(v(3,1),'%.2f') ' ' num2str(v(3,2),'%.2f') ' ' num2str(v(3,3),'%.2f') ']'],...
                ['L1/L2=' num2str(l(1)/l(2),'%1.1e\n') ', L2/L3=' num2str(l(2)/l(3),'%1.1e\n')],...  
                ['\theta_B=' num2str(ang2B,'%.0f') '^o, \theta_x=' num2str(ang2k,'%.0f') '^o']};
    title(hca,titlestr)
end
end
%% make time series to do minvar analysis
yind =13;
Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) tocolumn(squeeze(BZ(yind,:))*1e9)];
%By = [tocolumn(t) tocolumn(squeeze(BY(yind,:))*1e9)];
%Bz = [tocolumn(t) tocolumn(squeeze(BZ(yind,:))*1e9)];

%% check how angle changes when the y-location changes
for yind=1:numel(y);
    Bt = [tocolumn(t) tocolumn(squeeze(BX(yind,:))*1e9) tocolumn(squeeze(BY(yind,:))*1e9) tocolumn(squeeze(BZ(yind,:))*1e9)];
    [~,l,v]=irf_minvar(Bt);
    ang2B_vec(yind)=acosd(abs(([0 0 1]*v(3,:)')));
    ang2x_vec(yind)=acosd(abs(([1 0 0]*v(3,:)')));
end
if 1
%plot(ang2B_vec,y*1e-3)
%ylabel('y [km]')
figure(4)
npl=2;
for kk=1:npl; h(kk)=subplot(1,npl,kk); end
isub = 1;

hca=h(isub); isub=isub+1;
plot(hca,ang2B_vec,y*1e-3,ang2x_vec,y*1e-3)
xlabel(hca,'\theta  (acosd(|v_3 dot vector|)) [degrees]')
ylabel(hca,'y [km]')
irf_legend(hca,{'\theta_B','\theta_x'},[1 0.95])
set(hca,'ylim',y([1 end])*1e-3,'xlim',[0 90])

hca=h(isub); isub=isub+1;
hold(hca,'on')
tstop=fix(numel(t)/3);
ply = 1:3:tstop; % plot step
plx = 1:3:numel(y); % plot step
contour(hca,fX(:,1:tstop)*1e-3,fY(:,1:tstop)*1e-3,BZ(:,1:tstop))
quiver(hca,fX(plx,ply)*1e-3,fY(plx,ply)*1e-3,BX(plx,ply),BY(plx,ply))
ylabel(hca,'y [km]')
xlabel(hca,'x (=vph*t) [km]')
title(h(1),['LHDW B reconstruction, \theta = ' num2str(prop_angle*180/pi) '^o',...
            ', v_{ph}=' num2str(vph*1e-3) ' km/s'])
box(hca,'on')
%axis equal
set(hca,'ylim',y([1 end])*1e-3,'xlim',vph*t([1 tstop])*1e-3)

end
%ylabel('y [km]')
     %%
figure(3)
npl=5;
for kk=1:npl; h(kk)=subplot(npl,1,kk); end
isub = 1;

hca=h(isub); isub=isub+1;
plot(hca,t,phi(k,0,w,t))
ylabel(hca,'\phi [V]')

hca=h(isub); isub=isub+1;
plot(hca,t,E(kper,0,w,t,0)*1e3,t,E(kpar,0,w,t,0)*1e3)
ylabel(hca,'E [mV/m]')
irf_legend(hca,{'E_{\perp}','E_{||}'},[1 0.95])

hca=h(isub); isub=isub+1;
plot(hca,t,v_par(k,0,w,t,pi/2)*1e-6) % 10^3 km/s
ylabel(hca,'v_{e,||} [10^3 km/s]')

hca=h(isub); isub=isub+1;
plot(hca,t,-n*e*v_par(kpar,0,w,t,pi/2)*1e9); % nA
ylabel(hca,'j_{e,||} [nA]')

%hca=h(isub); isub=isub+1;
%plot(hca,t,-n*e*v_par(k,0,w,t,90)*1e9./phi(k,0,w,t)); % nA
%ylabel(hca,'j_{e,||}/\phi [nA/V]')


