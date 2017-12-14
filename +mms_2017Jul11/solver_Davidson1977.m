%% Observed wave phase velocity and wavelength
v_timing = 1.34e+03*[-0.85 -0.48  0.24];
v_amplitude = sqrt(sum(v_timing.^2));
l_obs = 100e3;
k_obs = 2*pi/l_obs;

%% Set up
units = irf_units;

% Plasma properties
% n = 0.06e6;
% ne = n;
% ni = n;
% Te = 500; Te_K = Te*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
% Ti = 7000; Ti_K = Ti*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
% B = 10e-9; 

% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

% Physical parameters
vti = sqrt(2*qe*Ti/mi);
vte = sqrt(2*qe*Te/me);
wpe = sqrt(ne*qe^2/(me*eps)); % Hz
wpi = sqrt(ni*qe^2/(mi*eps)); % Hz
wce = qe*B/me;
wci = qe*B/mi;
roe = vte/wce;
roi = vti/wci;
Ld = sqrt(vte^2)/wpe/sqrt(2)*1e-3;
wlh = sqrt(wce*wci);
PB = B^2/2/units.mu0;
Pi = Ti_K*units.kB*n;
Pe = Te_K*units.kB*n;
betai = Pi/PB;
betae = Pe/PB;
beta = betai + betae;


% Gradient length scales, m
% E = -50*1e-3; % E (together with Ln/LT) is used to get Ln and other parameters
% LnLT = 0.1; % Ln/LT, Ln generally smaller then LT
vE = -E/B;
vdi = -vE;
vde = -Te/Ti*vdi;

Ln = -(1+LnLT)*(Te_K*units.kB/units.e/B)./vde;
LT = Ln/LnLT;
LB = -2/(beta/Ln + betae/LT); % Davidson1975 eq. 13: (1/LB) = 0.5*beta*(1/Ln) - 0.5*betae*(1/LT)

vB = -1/LB*vte^2/2/wce;
vn = -1/Ln*vte^2/2/wce;
vT = -1/LT*vte^2/2/wce;    

disp('----------------------------------------------')
fprintf('B = %g nT \nn = %g cc \nTe = %g eV \nTi = %g eV \nbeta = %g \nbeta_i = %g \nbeta_e = %g\nvti = %g km/s\nroe = %.1f km\n\n',B*1e9,n*1e-6,Te,Ti,beta,betai,betae,vti*1e-3,roe*1e-3)
fprintf('E = %g mV/m \nLn/LT = %g \n\n',E*1e3,LnLT)
fprintf('Ln = %.0f km \nLB = %.0f km \nLT = %.0f km \n\n',Ln*1e-3,LB*1e-3,LT*1e-3)
fprintf('vE = %.0f km/s \nvde = %.0f km/s \nvB = %.0f km/s \nvn = %.0f km/s \nvT = %.0f km/s \n',vE*1e-3,vde*1e-3,vB*1e-3,vn*1e-3,vT*1e-3)
  
% Plot distributions
if 0
figure(52)
fe77_drift = @(v,Ln,LT,LB,vE,vte) (1/pi/vte^2)^(3/2)*exp(-(v-vE-vde).^2/vte.^2);
fe77_exp = @(v,Ln,LT,LB,vE,vte) (1/pi/vte^2)^(3/2)*exp(-(v-vE).^2/vte.^2).*(1-1*Ln^-1.*(v-vE)/wce.*(1-Ln/LT.*(3/2-v.^2/vte^2)));
fi77 = @(v,vti) (1/pi/vti^2)^(3/2)*exp(-v.^2/vti.^2);
vvec = linspace(-1*vte,2*vte,1000);
hh = plotyy(vvec*1e-3,[fe77_exp(vvec,Ln,LT,LB,vE,vte)',fe77_drift(vvec,Ln,LT,LB,vE,vte)'],vvec*1e-3,fi77(vvec,vti));
hh(2).YLim = [0 max(fi77(vvec,vti))*1.1];

legend('expanded f_e','drifting maxwellian; ve = vE + vde',' ions stationary maxwellian')
end
%% Dispersion solver
% Wavenumber vector
%nk = 200;
% kmin= 1e-1; kmax = 2;
kvec = linspace(k_min,k_max,nk)/roe;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
x = 0;
for ik = 1:nk
  if guess_previous_k 
    xguess = x;
  else
    xguess = wr_all(iE,iTe-1,ik) + wi_all(iE,iTe-1,ik)*1i;
  end
 % af = @(temp) mms_2017Jul11.D_Davidson1977(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'colde');  
  af = @(temp) mms_2017Jul11.D_Davidson1977(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full',0);  
 % af = @(temp) mms_2017Jul11.D_Davidson1975(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);  
  options = optimset('GradObj','on','display','off');  
  [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
  %fprintf('x = %g + %gi   ',real(x),imag(x))
  %mms_2017Jul11.D_Davidson1977(x,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full',1);
  wr_store(ik) = real(x);
  wi_store(ik) = imag(x);  
end
 
ikmax = find(wi_store==max(wi_store));
kmax = kvec(ikmax);
vphmax = wr_store(ikmax)/kvec(ikmax);
wimax = wi_store(ikmax);
wrmax = wr_store(ikmax);

%% Plot solution
if 0
figure(103)
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

isub = 1;
if 1
  hca = h(isub); isub = isub + 1;  
  s1 = sprintf('B = %g nT \nn = %g cc \nTe = %g eV \nTi = %g eV \nbeta = %g \nbeta_i = %g \nbeta_e = %g\nvti = %g km/s\nroe = %.1f km',B*1e9,n*1e-6,Te,Ti,beta,betai,betae,vti*1e-3,roe*1e-3);
  s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \nvn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \nE = %.0f mV/m \n',Ln*1e-3,LB*1e-3,LT*1e-3,vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3,E*1e3);
  s3 = sprintf('k_{max}*roe = %g \nwi_{max}/wlh = %.3f \nwr_{max}/wlh = %.3f \nvph_{max} = %.0f km/s',kmax*roe,wimax/wlh,wrmax/wlh,vphmax*1e-3);

  if exist('ht1','var'); delete(ht1); end
  if exist('ht2','var'); delete(ht2); end
  if exist('ht3','var'); delete(ht3); end
  ht1 = text(hca,hca.XLim(1),hca.YLim(2),s1,'verticalalignment','top','horizontalalignment','left');
  ht2 = text(hca,mean(hca.XLim)*0.7,hca.YLim(2),s2,'verticalalignment','top','horizontalalignment','left');
  ht3 = text(hca,hca.XLim(2),hca.YLim(2),s3,'verticalalignment','top','horizontalalignment','right');
  hca.Visible = 'off';
end
if 1
  hca = h(isub); isub = isub + 1;
  ax = plotyy(hca,kvec*roe,wr_store/wlh,kvec*roe,wi_store/wlh);
  ax(1).XLabel.String = 'k\rho_e';
  ax(1).YLabel.String = '\omega/\omega_{LH}';
  ax(2).YLabel.String = '\gamma/\omega_{LH}';
   
  %ax(1).YLim = [0 max(wr_store)/wlh*1.2];
  %ax(2).YLim = [0 wimax/wlh*1.2];
  
  ax(2).YGrid = 'on';
  ax(2).YMinorGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,wr_store./kvec*1e-3,kmax*roe,wr_store(ikmax)/kmax*1e-3,'*k');
  hca.XLabel.String = 'k\rho_e';
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.YGrid = 'on';  
end
drawnow
end
%%
if 0
dispcurv = struct;
dispcurv.kvec = kvec;
dispcurv.wr = wr_;
dispcurv.wi = wi_; 
dispcurv.wpi = wpi; 
dispcurv.wpe = wpe;
dispcurv.R = 0.25;
dispcurv.ne1 = ni*(1-R);
dispcurv.ne2 = ni*R;
dispcurv.Te1 = 300;
dispcurv.Te2 = 12;
dispcurv.Ti = 300;
dispcurv.S = 0.5;
dispcurv.Ld = Ld;
dispcurv.veth1 = veth1;
dispcurv.veth2 = veth2;
dispcurv.vith = vith;

if 0,
save('modbun.mat','dispcurv');
end

end