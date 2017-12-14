% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

%% Set up
units = irf_units;

% Plasma properties
n = 0.07e6;
ne = n;
ni = n;
Te = 500; Te_K = Te*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti = 5000; Ti_K = Ti*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
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

% gradient drifts
vB = 1/LB*vte^2/2/wce;
vn = -1/Ln*vte^2/2/wce;
vT = -1/LT*vte^2/2/wce;

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

%% Examine dispersion relation
% besseli Modified Bessel function of the first kind
if 0
J0p = @(x) (s*besselj(0, x))./x - besselj(0 + 1, x); % x is mu = k*v/wce
J0 = @(x) besselj(0, x); % x is mu = k*v/wce
A = @(w,k,v) w - k*vE - k*vn - k*vT*(1-v.^2/vte^2);

O1_ = @(w,k,v) 2/vte^4*v.^3.*J0p(k*v/wce).^2          .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
O2_ = @(w,k,v) 2/vte^3*v.^2.*J0(k*v/wce).*J0p(k*v/wce).*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
O3_ = @(w,k,v) 2/vte^2*v.^1.*J0(k*v/wce).^2           .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);

nv = 10000;
vvec = linspace(1,5*vte,nv);

ik = 50;
k = kvec(ik);
w = wr_(ik) + i*wi_(ik);

O1 = trapz(O1_(wtot,k,vvec),vvec);
O2 = trapz(O2_(wtot,k,vvec),vvec);
O3 = trapz(O3_(wtot,k,vvec),vvec);
O1, O2, O2

D = mms_2017Jul11.D_Davidson1977(w,k,vte,wce,wpe,vti,wci,wpi,vE,vB,vn,vT);
end

%%
vE = vE-0*vB;
  
% Observed wave phase velocity and wavelength
v_timing = 1.34e+03*[-0.85 -0.48  0.24];
v_amplitude = sqrt(sum(v_timing.^2));
l_obs = 100e3;
k_obs = 2*pi/l_obs;
roek_obs = k_obs*roe;



%mms_2017Jul11.LH_Davidson1975_plot_dispersion_properties;

% Dispersion solver
nk = 50;
kvec = 1*linspace(1e-2,2,nk)/roe;
dk = kvec(2)-kvec(1);
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));


%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
% higher vde2 -> higher vph, this should control the guessing range
% w/k = vph -> wrange >~vph*k

vguess = 1*vE;
wapprox = 1*kvec(1)*vguess;%wlh;kvec(1)*vguess;
wrvec = linspace(0,1,1501)*wapprox;
wivec = linspace(0,1,1501)*wapprox;
[wr,wi] = meshgrid(wrvec,wivec);
wmat = wr+i*wi;

do75 = 0;
D75 = mms_2017Jul11.D_Davidson1975(wmat,kvec(1),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
D77 = mms_2017Jul11.D_Davidson1977(wmat,kvec(1),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
if do75
  D_ = D75;
else
  D_ = D77;
end
z = abs(D_); 
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

figure(17)

h(1) = subplot(3,1,1);
h(2) = subplot(3,1,2);
h(3) = subplot(3,1,3); h(3).Box='on';
plot(h(1),kvec(1)*roe,wr_(1)/wlh,'o',kvec(1)*roe,wi_(1)/wlh,'x'); 


colors = mms_colors('matlab');
if exist('hl','var'); delete(hl); end
hl = plot(h(2),kvec([1 end])*roe,v_amplitude*[1 1]...,'color',colors(1,:)...
              ,k_obs*roe,v_amplitude,'*'...%,'color',colors(2,:)...
              ...%,kvec*roe,v_delta(kvec)*1e-3...%,'color',colors(3,:)...
              ,kvec([1 end])*roe,vE*1e-3*[1 1]...%,'color',colors(4,:)...
          );        

hl(1).Color = colors(1,:);
hl(2).Color = colors(1,:);
hl(3).Color = colors(2,:);
%hl(4).Color = colors(3,:);

legend(hl,'vph obs','k obs','vE','v\Delta')        

colors = mms_colors('matlab');

x_real_store = [];
x_imag_store = [];
x_real_store(1) = 0;
x_imag_store(1) = 0;

frange = 1.0*wapprox;
for nn = 2:length(kvec)
  frange = vguess*dk*2;
  %linear search for zero point around linearly interpolated guess. frange
  %may need to be changed depending on mode. Large frange decreases accuracy
  nf = 401;
%   wrv = linspace(wr_(nn-1)+gradr-frange,wr_(nn-1)+gradr+frange,nf+1);
%   wiv = linspace(wi_(nn-1)+gradi-frange,wi_(nn-1)+gradi+frange,nf);
  wrv = linspace(0,5*wlh,nf);
  wiv = linspace(0,2*wlh,nf);
  %if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end
  isneg = find(wrv<0);
  wrv(isneg)=[];
  %wiv(isneg)=[];

  [wr,wi] = meshgrid(wrv,wiv);
  wc = wr+i*wi; 
  wc = wc;
  
  if 1 % compare Davidson1975 and Davidson1977
    D75 = mms_2017Jul11.D_Davidson1975(wc,kvec(nn),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
    D77 = mms_2017Jul11.D_Davidson1977(wc,kvec(nn),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
    %D77 = D75;
    z75 = abs(D75);
    z77 = abs(D77);
    
    figure(13)
    nrows = 3;    
    ncols = 2;
    npanels = nrows*ncols;
    for ip = 1:npanels
      h2(ip) = subplot(nrows,ncols,ip);
    end
    isub = 1;
    if 1 % Davdison1975     
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,abs(D75)); shading(hca,'flat');
      hca.Title.String = '1975, abs(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(z75)) max(max(z75))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,abs(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, abs(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(z75)) max(max(z75))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1975     
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,real(D75)); shading(hca,'flat');
      hca.Title.String = '1975, real(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(real(D75))) max(max(real(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,real(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, real(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(real(D77))) max(max(real(D77)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1975     
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,imag(D75)); shading(hca,'flat');
      hca.Title.String = '1975, imag(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(imag(D75))) max(max(imag(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrv/wlh,wiv/wlh,imag(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, imag(D)';
      colorbar('peer',hca)
      hca.CLim = [min(min(imag(D75))) max(max(imag(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    %pause
  end
  if do75
    z = z75;
  else
    z = z77;
  end  
  
  residue = min(min(z));
  [q r] = find(min(min(z))==z);
  wr_(nn) = wrv(r(1));
  wi_(nn) = wiv(q(1));
  
  
  af = @(temp) mms_2017Jul11.D_Davidson1977(temp,kvec(nn),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
  %[x,FVAL,EXITFLAG] = fsolve(af,[x_real_store(end) x_imag_store(end)], optimset('GradObj','on','display','off'));
  %[x,FVAL,EXITFLAG] = fsolve(af,[wr_(nn)  1*wi_(nn)], optimset('GradObj','on','display','off'));
  %[x,FVAL,EXITFLAG] = fsolve(af,wr_(nn)+i*wi_(nn), optimset('GradObj','on','display','off'));
  [x,FVAL,EXITFLAG] = fsolve(af,wlh*[1+1*i]*linspace(0,1,10), optimset('GradObj','on','display','off'));
  x
  FVAL
  minFVAL = min(abs(FVAL));
  iminFVAL = find(abs(FVAL)==minFVAL);
  x_real_store = [x_real_store real(x(iminFVAL))]; 
  x_imag_store = [x_imag_store imag(x(iminFVAL))];
  
  wnorm = wlh;
  if 1
    figure(17)
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
    plot(hca,kvec(nn)*roe,real(x)/wnorm,'+b',kvec(nn)*roe,imag(x)/wnorm,'+r')
    plot(hca,kvec(nn)*roe,x_real_store(nn)/wnorm,'+b',kvec(nn)*roe,x_imag_store(nn)/wnorm,'+r')
    hca.XLabel.String = 'k*\rho_e';
    hca.YLabel.String = '\gamma/\omega_{LH}, \omega/\omega_{LH}';
    hca.YLim = [-2 2];
    %hca.YScale = 'log';
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
  pause%(0.1)
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