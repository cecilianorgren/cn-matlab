% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

units = irf_units;

% Plasma properties
n = 0.07e6;
ne = n;
ni = n;
Te = 500; Te_K = Te*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti = 4000; Ti_K = Ti*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
B = 10e-9;
 
% new better ion temperature using gsePi and ne
%c_eval('gseTi?_ = irf.ts_tensor_xyz(gsePi?.time,gsePi?.data*1e-9./repmat(ne?.resample(gsePi?).data*1e6,1,3,3)/units.kB*units.kB/units.e);',1)
%figure(22)
%irf_plot({gseTi1.trace/3,gseTi1_.trace/3},'comp')

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
Ln = 200e3;
LT = 800e3;
LB = 300e3; % can be overridden below if self-consistent approach is chosen
self_consistent_LB = 1;

% Davidson1975 eq. 13
% (1/LB) = 0.5*beta*(1/Ln) - 0.5*betae*(1/LT)
if self_consistent_LB
  LB = 2/(beta/Ln + betae/LT);
end

% Cross field drifts
vE = 2000e3; % can be overridden below if self-consistent approach is chosen
self_consistent_vE = 1;
vdi = 0;

plot_v_delta = 0;

% Consistency check

disp('----------------------------------------------')
% input parameters
fprintf('B = %g nT \nn = %g cc \nTe = %g eV \nTi = %g eV \nbeta = %g \nbeta_i = %g \nbeta_e = %g\nvti = %g km/s\n\n',B*1e9,n*1e-6,Te,Ti,beta,betai,betae,vti*1e-3)
fprintf('Ln = %g km \nLB = %g km \nLT = %g km \n\n',Ln*1e-3,LB*1e-3,LT*1e-3)

% magnetic field gradient drift
vB = 1/LB*vte^2/2/wce;

disp('vB from LB:')
fprintf('vB = vte^2/(LB*2*wce) = %.0f km/s \n',vB*1e-3)

% vde from LN, LT, B
vde_L = (Te_K*units.kB/units.e/B)*(1/Ln + 1/LT);
disp(' ')
disp('vde from Ln and LT:')
fprintf('vde = Te*kB/eB*(1/Ln+1/LT) = %.0f km/s \n',vde_L*1e-3)

disp(' ')
if self_consistent_vE
  vE = vde_L*Ti/Te;
  disp('vE from gradient length scales and ion force balance:')
  fprintf('vE = -vdi = Ti/Te*vde = Te*kB/eB*(1/Ln+1/LT)*Ti/Te = %g km/s \n',vE*1e-3)
else
  disp('vE as free paramater:')
  fprintf('vE (free parameter) = %g km/s \n',vE*1e-3)
end

% electric field corresponding to vE
E = vE*B;
disp('E corresponding to vE:')
fprintf('E = vE*B = %g mV/m \n',E*1e3)
  
% for force balance on fluid element, vExB = -vdi = (Ti/Te)*vde
% ve = vExB + Vde = vExB - (Te/Ti)*vdi = (1+Te/Ti)*vExB
% vde = (Te/Ti)*vExB
vde_E = (Te/Ti)*vE;

disp(' ')
disp('vde from ion force balance:')
fprintf('vde = -(Te/Ti)*vdi = (Te/Ti)*vExB = %.0f km/s \n',vde_E*1e-3)

disp(' ')
if self_consistent_LB
  disp('LB self consistent: LB = (0.5*beta*(1/Ln) - 0.5*betae*(1/LT))^-1')
else
  disp('LB as free parameter:')
end
disp('ratio of gradient length scales: (Davidson 1975, (1/LB) = 0.5*beta*(1/Ln) - 0.5*betae*(1/LT))')  
fprintf('0.5*beta*LB/Ln = %g \n',0.5*beta*LB/Ln)
fprintf('0.5*betae*LB/LT = %g \n',0.5*betae*LB/LT)


vE = vE-0*vB;
  
%
% Observed wave phase velocity and wavelength
v_timing = 1.34e+03*[-0.85 -0.48  0.24];
v_amplitude = sqrt(sum(v_timing.^2));
l_obs = 100e3;
k_obs = 2*pi/l_obs;
roek_obs = k_obs*roe;

% Dispersion relation
Nfad = 50;
xi = @(w,k,vd,vt) (w-k*vd)/(k*vt);
b = @(k,vt,wc) k.^2*vt^2/2/wc^2;

% v_delta: Davidson1975, eq. 7
if 1
  % only Ln
  v_delta_Ln = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...
        +1/Ln);
  % Ln and LT
  v_delta_LnLT = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...
        +1/Ln...
        -1/LT*b(k,vte,wce).*(1-besseli(1,b(k,vte,wce))./besseli(0,b(k,vte,wce))));
  % Ln LT and LB
  v_delta_LnLBLT = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...
        +1/Ln...
        -1/LB*(1-b(k,vte,wce).*(1-besseli(1,b(k,vte,wce))./besseli(0,b(k,vte,wce))))...
        -1/LT*b(k,vte,wce).*(1-besseli(1,b(k,vte,wce))./besseli(0,b(k,vte,wce))));
   % only LB   
   v_delta_LB = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...        
        -1/LB*(1-b(k,vte,wce).*(1-besseli(1,b(k,vte,wce))./besseli(0,b(k,vte,wce)))));   
  % only LT
  v_delta_LT = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...        
        -1/LT*b(k,vte,wce).*(1-besseli(1,b(k,vte,wce))./besseli(0,b(k,vte,wce))));
  
  if plot_v_delta
    nk = 80; % just for plotting v_delta, changed below for the dispersion solver
    kvec = 1*linspace(1e-3,1,nk)/roe; 
  
    figure(77)    
    hca = subplot(1,1,1);    
    plot(hca,...
         kvec*roe,v_delta_Ln(kvec)*1e-3,...
         kvec*roe,v_delta_LT(kvec)*1e-3,...
         kvec*roe,v_delta_LB(kvec)*1e-3,...
         kvec*roe,v_delta_LnLT(kvec)*1e-3,...
         kvec*roe,v_delta_LnLBLT(kvec)*1e-3)
    legend(hca,'Ln','LT','LB','Ln,LT','Ln,LB,LT')
    hca.XLabel.String = 'k\rho_e';
    hca.YLabel.String = 'v_\Delta (km/s)';
  end
  v_delta = v_delta_LnLBLT;
end


%psi_i = @(w,k,vd,vt,wp) 2*wp^2/(k^2*vt^2)*(1 + i*sqrt(pi)*xi(w,k,vd,vt).*faddeeva(xi(w,k,vd,vt),Nfad));
psi_i = @(w,k,vt,wp,wc) 2*wp^2./(k.^2)/(vti^2).*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad));
psi_e = @(w,k,vt,wp,wc) (wp^2/wc^2).*(1-besseli(0,b(k,vt,wc)).*exp(-b(k,vt,wc)))./b(k,vt,wc) + ...
                           2*wp^2./k.^2/vt^2.*k.*v_delta(k)./(w-k*vE);

D = @(w,k) 1 + psi_i(w,k,vti,wpi) + psi_e(w,k,vte,wpe,wce);

%mms_2017Jul11.LH_Davidson1975_plot_dispersion_properties;

% Dispersion solver
nk = 100;
kvec = 0.5*linspace(1e-2,1,nk)/roe;
dk = kvec(2)-kvec(1);
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));


%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
% higher vde2 -> higher vph, this should control the guessing range
% w/k = vph -> wrange >~vph*k

vguess = 1*vE;
wapprox = kvec(1)*vguess;wlh;kvec(1)*vguess;
wrvec = linspace(0,1,1501)*wapprox;
wivec = linspace(0,1,1501)*wapprox;
[wr,wi] = meshgrid(wrvec,wivec);
wmat = wr+i*wi;

z = abs(D(wmat,kvec(1))); 
%min(min(z));
%pcolor(abs(z)); shading flat
[q r] = find(min(min(z))==z);
wr_(1) = wrvec(r);
wi_(1) = wivec(q);
 
gradr = wr_(1);
gradi = wi_(1);

if 0 % check grid    
    %%
    figure(79)
    subplot(3,1,1)
    pcolor(wr,wi,real(D(wmat,kvec(1))))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,2)
    pcolor(wr,wi,imag(D(wmat,kvec(1))))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,3)
    pcolor(wr,wi,abs(D(wmat,kvec(1))))
    xlabel('w_r')
    ylabel('w_i')
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[0 1])
end

figure(78)

h(1) = subplot(3,1,1);
h(2) = subplot(3,1,2);
h(3) = subplot(3,1,3); h(3).Box='on';
plot(h(1),kvec(1)*roe,wr_(1)/wlh,'o',kvec(1)*roe,wi_(1)/wlh,'x'); 


colors = mms_colors('matlab');
if exist('hl','var'); delete(hl); end
hl = plot(h(2),kvec([1 end])*roe,v_amplitude*[1 1]...,'color',colors(1,:)...
              ,k_obs*roe,v_amplitude,'*'...%,'color',colors(2,:)...
              ,kvec*roe,v_delta(kvec)*1e-3...%,'color',colors(3,:)...
              ,kvec([1 end])*roe,vE*1e-3*[1 1]...%,'color',colors(4,:)...
          );        

hl(1).Color = colors(1,:);
hl(2).Color = colors(1,:);
hl(3).Color = colors(2,:);
hl(4).Color = colors(3,:);

legend(hl,'vph obs','k obs','v\Delta','vE')        

colors = mms_colors('matlab');
frange = 1.0*wapprox;
for nn = 2:length(kvec)
  frange = vguess*dk*2;
  %linear search for zero point around linearly interpolated guess. frange
  %may need to be changed depending on mode. Large frange decreases accuracy
  nf = 201;
  wrv = linspace(wr_(nn-1)+gradr-frange,wr_(nn-1)+gradr+frange,nf);
  wiv = linspace(wi_(nn-1)+gradi-frange,wi_(nn-1)+gradi+frange,nf);
  %if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end
  %wrv(wrv<0)=[];

  [wr,wi] = meshgrid(wrv,wiv);
  wc = wr+i*wi;

  z = abs(D(wc,kvec(nn)));
  ze = abs(psi_e(wc,kvec(nn),vte,wpe,wce));
  zi = abs(psi_i(wc,kvec(nn),vti,wpi,wci));
  %z=ze;
  residue = min(min(z));
  [q r] = find(min(min(z))==z);
  wr_(nn) = wrv(r);
  wi_(nn) = wiv(q);

  wnorm = wlh; wpi;
  if 1
    hca = h(1);
    hold(hca,'on')
    if min(min(z))>1e-3
      color1 = colors(3,:);
      color2 = colors(3,:);
    else
      color1 = colors(1,:);
      color2 = colors(2,:);
    end
    plot(hca,kvec(nn)*roe,wr_(nn)/wnorm,'.b',kvec(nn)*roe,wi_(nn)/wnorm,'.r')
    hca.XLabel.String = 'k*\rho_e';
    hca.YLabel.String = '\gamma/\omega_{LH}, \omega/\omega_{LH}';
    %errorbar(h(1),kvec(nn)*roe,wrv(fix(nf/2))/wnorm,frange/wnorm)
    %errorbar(h(1),kvec(nn)*roe,wiv(fix(nf/2))/wnorm,frange/wnorm)
    
    hca = h(2);
    hold(hca,'on')
    plot(hca,kvec(nn)*roe,wr_(nn)/kvec(nn)*1e-3,'.k')
    hca.XLabel.String = 'k*\rho_e';
    hca.YLabel.String = 'v_{ph} (km/s)';
    
    hca = h(3);
    hold(hca,'on')
    semilogy(hca,kvec(nn)*roe,residue,'.')
    hca.XLabel.String = 'k*\rho_e';
    hca.YLabel.String = 'residue';
    hca.YScale = 'log';
  else
    pcolor(abs(z)); shading flat;
    colorbar
    set(gca,'clim',[0 1e2])
  end
  pause(0.1)
  %To do: Write better root finder.

  gradr = wr_(nn)-wr_(nn-1);
  gradi = wi_(nn)-wi_(nn-1);
end

%%
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

pos=find(max(wga)==wga);
kmax = kvec(pos(1));
vph = wia(pos(1))/kvec(pos(1));
vphnorm = vph(1)/veth1;
maxfreq = wia(pos(1))/wpi;
maxgamma = wga(pos(1))/wpi;

if 0
%%

plot(kvec*Ld,wr_/wpi,kvec*Ld,wi_/wpi)
plot(kvec*Ld,wr_/wpe,kvec*Ld,wi_/wpe)
xlabel('k\lambda_D')
ylabel('\omega_{max}/\omega_{pi}, \gamma_{max}/\omega_{pi}')
distrStr = {['R=' num2str(ne2/ni,'%.2f')],...
            ['Ti=' num2str(Ti,'%.0f')],...
            ['Te1=' num2str(Te1,'%.0f')],...
            ['Te2=' num2str(Te2,'%.0f')],...
            ['S=' num2str(vd2/veth1,'%.2f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')
end