
mu0 = 1.2566e-06;
M = 8.22e22; % A m^2
RE = 6.3712e+03; % km
h = 100;

r0 = RE;
th0 = 90-60; % at field line crossing Earth's surface, from lat to colat
r1 = RE + h;
th2 = 90; % at equatorial plan

% We have r1 = Re+h, th2=90 deg, we need to get th1, r2 

th1 = asind(sqrt((RE+h)/4/RE));
%th1 = 30; % trying to see if assuming th2 = th0 changes anything?
r2 = r0*sind(th2)^2/sind(th0)^2;


f_B = @(r,th) mu0*M/(4*pi*r^3)*sqrt(1+3*cosd(th)^2);
f_Bth = @(th) mu0*M/(4*pi*(4*RE)^3)*sqrt(1+3*cosd(th)^2)/sind(th)^6;

B1 = f_B(r1,th1);
B2 = f_B(r2,th2);

pa1 = 90; % mirroring

pa2 = asind(sqrt(B2/B1)*sind(pa1))

B2B1 = (sqrt(1+3*cosd(th2))/sind(th2)^6)/(sqrt(1+3*cosd(th1))/sind(th1)^6);
B2B1_th90 = (sqrt(1+3*0)/1^6)/(sqrt(1+3*cosd(th1))/sind(th1)^6);

% c
h = 100; % km
th0 = 90-[15 30 45 60 70];
th1 = asind(sqrt((RE+h)/RE)*sind(th0));

B2B1 = sind(th1).^6./sqrt(1+3*cosd(th1).^2);
pa2 = asind(sqrt(B2B1));

%% Plasma pause radius
mu0 = 1.2566e-06;
RE = 6.3712e+06; % m
M = 8.22e22; % A m^2
wE = 2*pi/(24*60*60); % Earth's rotation frequency, rad/s
Econv = 1*[0.5 2]*1e3/(RE); % kV/RE - > V/m (to have V in SI units)
r = (mu0*M*wE./(4*pi*Econv)).^(1/2);

Bz = @(r) mu0*M./(4*pi*r.^3);
vcorot = @(r) wE.*r;
Ecorot = @(r) Bz(r).*vcorot(r);
vconv = @(r,Econv) Econv./Bz(r);

rmin = 0; rmax = 15*RE;
rvec = linspace(rmin,rmax,100);

colors = [     0    0.4470    0.7410;...
          0.8500    0.3250    0.0980;...
          0.9290    0.6940    0.1250];
clear h
nrows = 3; ncols = 1; ipanel = 0;
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

hca = h(isub); isub = isub + 1;
semilogy(hca,rvec/RE,Bz(rvec))
hca.YTick = 10.^(-7:1:10);
hca.YLabel.String = 'Magnetic field (T)';
hca.XLabel.String = 'Radius (R_E)';

hca = h(isub); isub = isub + 1;
semilogy(hca,rvec/RE,vcorot(rvec))
hca.YTick = 10.^(-7:1:10);
hca.YLabel.String = 'Speed (m/s)';
hca.XLabel.String = 'Radius (R_E)';

hca = h(isub); isub = isub + 1;
semilogy(hca,rvec/RE,Ecorot(rvec),...
  [rmin rmax]/RE,Econv(1)*[1 1],...
  [rmin rmax]/RE,Econv(2)*[1 1])
hold(hca,'on')
set(hca,'ColorOrder',colors(2:3,:))
semilogy(r(2)/RE,Econv(2),'*')
semilogy(r(1)/RE,Econv(1),'*')
hold(hca,'off')
hca.YTick = 10.^(-7:1:10);
hca.YLabel.String = 'Electric field (V/m)';
hca.XLabel.String = 'Radius (R_E)';

c_eval('h(?).FontSize = 14;',1:numel(h))
c_eval('h(?).XTick = [0:2:20];',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
compact_panels(0.02,0.00)