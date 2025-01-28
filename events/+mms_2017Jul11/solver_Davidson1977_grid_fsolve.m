% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

%% Observed wave phase velocity and wavelength
v_timing = 1.34e+03*[-0.85 -0.48  0.24];
v_amplitude = sqrt(sum(v_timing.^2));
l_obs = 100e3;
k_obs = 2*pi/l_obs;

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

input = 'E';
switch input % 'Ln', 'E';
  case 'Ln'
    %%
    % Gradient length scales, m
    Ln = 450e3;
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
    vn = 1/Ln*vte^2/2/wce;
    vT = 1/LT*vte^2/2/wce;

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
  case 'E'  
    %%
    % Gradient length scales, m
    E = -10*1e-3; % V/m
    LnLT = 0.1; % Ln/LT, Ln generally smaller then LT
    
    vE = -E/B;
    vdi = -vE;
    vde = -Te/Ti*vdi;
    
    Ln = -(1+LnLT)*(Te_K*units.kB/units.e/B)/vde;
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
end

%% Set up figure
doplot = 1;
if doplot
  colors = mms_colors('matlab');
  figure(17)
  nrows = 4;    
  ncols = 1;
  npanels = nrows*ncols;
  for ip = 1:npanels
    h(ip) = subplot(nrows,ncols,ip);
  end

  if 0
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

  legend(hl,'vph obs','k obs','vE','v\Delta')  
  end
  %linkaxes(h,'x')
end

%% Dispersion solver
% Wavenumber vector
nk = 100;
kvec = linspace(1e-1,2.5,nk)/roe;

% save vector
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));
resa = nan(1,length(kvec));

% Two step solver. 
% 1. Search a large interval of wr and wi to find approximate minimum.
% 2a. Use fsolve to get more precise solution
% 2b. Make smaller finer grid to get more precise solution

nw = 501;
wirange = 1*wlh;
wrrange = 3*wlh;
wrcenter = wrrange/2;
wicenter = wirange/2;

x_real_store = nan(1,nk);
x_imag_store = nan(1,nk);

for ik = 1:length(kvec)  
  % Initial grid of wr wi to find minimum of D
  wrvec = linspace(wrcenter-wrrange/2,wrcenter+wrrange/2,nw);
  wivec = linspace(wicenter-wirange/2,wicenter+wirange/2,nw);
  [wr,wi] = meshgrid(wrvec,wivec);
  wmat = wr + i*wi;

  %do75 = 0;
  %D75 = mms_2017Jul11.D_Davidson1975(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
  if ik > 1
    x_tmp = x;
  else
    tic
  D77 = mms_2017Jul11.D_Davidson1977(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full');
  toc
  %if do75, z = abs(D75);
  %else,    z = abs(D77);
  %end
  z = abs(D77);
  residue = min(min(z));
  [q r] = find(z==residue);
  wr_tmp(ik) = wrvec(r);
  wi_tmp(ik) = wivec(q);
  
  x_tmp = [wrvec(r) + 1i*wivec(q)];
  end
  if 0*doplot % compare Davidson1975 and Davidson1977
    %D75 = mms_2017Jul11.D_Davidson1975(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'colde');
    %D75 = mms_2017Jul11.D_Davidson1977(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'colde');
    %D77 = mms_2017Jul11.D_Davidson1977(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full');
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
      pcolor(hca,wrvec/wlh,wivec/wlh,abs(D75)); shading(hca,'flat');
      hca.Title.String = '1975, abs(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(z75)) max(max(z75))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrvec/wlh,wivec/wlh,abs(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, abs(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(z75)) max(max(z75))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1975     
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrvec/wlh,wivec/wlh,real(D75)); shading(hca,'flat');
      hca.Title.String = '1975, real(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(real(D75))) max(max(real(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrvec/wlh,wivec/wlh,real(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, real(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(real(D77))) max(max(real(D77)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1975     
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrvec/wlh,wivec/wlh,imag(D75)); shading(hca,'flat');
      hca.Title.String = '1975, imag(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(imag(D75))) max(max(imag(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    if 1 % Davdison1977, finite beta
      hca = h2(isub); isub = isub + 1;
      pcolor(hca,wrvec/wlh,wivec/wlh,imag(D77)); shading(hca,'flat');
      hca.Title.String = '1977, finite beta, imag(D)';
      colorbar('peer',hca)
      %hca.CLim = [min(min(imag(D75))) max(max(imag(D75)))];
      hca.YLabel.String = 'w_i/wlh';
      hca.XLabel.String = 'w_r/wlh';
    end
    pause
  end
  
  %af = @(temp) mms_2017Jul11.D_Davidson1977(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full');
  af = @(temp) mms_2017Jul11.D_Davidson1977(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full');
  %af = @(temp) mms_2017Jul11.D_Davidson1975(temp,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT);
  options = optimset('GradObj','on','display','off'); 
 
  [x,FVAL,EXITFLAG] = fsolve(af,[wr_tmp(ik) + 1i*wi_tmp(ik)]*1,options);  
  
  disp(sprintf('%g + %gi, fval = %g + %gi',real(x),imag(x),real(FVAL),imag(FVAL)))

  minFVAL = min(abs(FVAL));
  iminFVAL = find(abs(FVAL)==minFVAL);
  resa(ik) = minFVAL;
  x_real_store(ik) = real(x(iminFVAL));
  x_imag_store(ik) = imag(x(iminFVAL));
  
  wnorm = wlh;
  if doplot % plot
    figure(17)
    isub = 1;    
    if 1 % wr
      hca = h(isub); isub = isub + 1;
      hold(hca,'on')      
      plot(hca,kvec(ik)*roe,real(x)/wnorm,'+r',kvec(ik)*roe,wr_tmp(ik)/wnorm,'.b')
      hca.XLabel.String = 'k*\rho_e';
      hca.YLabel.String = '\omega/\omega_{LH}';
    end
    if 1 % wi
      hca = h(isub); isub = isub + 1;
      hold(hca,'on')      
      plot(hca,kvec(ik)*roe,imag(x)/wnorm,'+r',kvec(ik)*roe,wi_tmp(ik)/wnorm,'.b')
      hca.XLabel.String = 'k*\rho_e';
      hca.YLabel.String = '\gamma/\omega_{LH}';
      
    end
    if 1 % wi
      hca = h(isub); isub = isub + 1;
      hold(hca,'on')
      plot(hca,kvec(ik)*roe,real(x)/kvec(ik)*1e-3,'+r',kvec(ik)*roe,wr_tmp(ik)/kvec(ik)*1e-3,'.k')
      hca.XLabel.String = 'k*\rho_e';
      hca.YLabel.String = 'v_{ph} (km/s)';
    end
    if 1 % residue abs(D)
      hca = h(isub); isub = isub + 1;
      hold(hca,'on')
      semilogy(hca,kvec(ik)*roe,residue,'.b',kvec(ik)*roe,abs(mms_2017Jul11.D_Davidson1977(x,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full')),'+r')
      hca.XLabel.String = 'k*\rho_e';
      hca.YLabel.String = 'residue';
      hca.YScale = 'log';
    end
  elseif 0
    pcolor(abs(z)); shading flat;
    colorbar
    set(gca,'clim',[0 1e2])
  end
  pause(0.1)

  if ik == 1    
    gradr = wr_tmp(1);
    gradi = wi_tmp(1);
  else
    gradr = wr_tmp(ik)-wr_tmp(ik-1);
    gradi = wi_tmp(ik)-wi_tmp(ik-1);
  end
  if ik == 1    
    gradr = x_real_store(1);
    gradi = x_imag_store(1);
  else
    gradr = x_real_store(ik)-x_real_store(ik-1);
    gradi = x_imag_store(ik)-x_imag_store(ik-1);
  end
 
  wrcenter = wrcenter + gradr;
  wicenter = wicenter + gradi;
  %wrcenter = real(x);
  %wicenter = imag(x);
  if 0%gradr > 20 % jump in solution, use finer search grid instead of fsolve
    wrcenter = x_real_store(ik-1);
    wicenter = x_imag_store(ik-1);
    wrrange = wlh/4;
    wirange = wlh/4;
    wrvec = linspace(wrcenter-wrrange/2,wrcenter+wrrange/2,nw);
    wivec = linspace(wicenter-wirange/2,wicenter+wirange/2,nw);
    [wr,wi] = meshgrid(wrvec,wivec);
    wmat = wr + i*wi;
        
    D77 = mms_2017Jul11.D_Davidson1977(wmat,kvec(ik),vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,'full');
    if do75, z = abs(D75);
    else,    z = abs(D77);
    end

    residue = min(min(z));
    [q r] = find(z==residue);
    wr_tmp(ik) = wrvec(r);
    wi_tmp(ik) = wivec(q);

    x_tmp =  [wrvec(r) + 1i*wivec(q)];
  end
end

%%
figure(103)
ikmax = find(abs(x_imag_store)==max(abs(x_imag_store)));
kmax = kvec(ikmax);
vph = x_real_store(ikmax)/kvec(ikmax);
wimax = x_imag_store(ikmax);
wrmax = x_real_store(ikmax);


h = plotyy(kvec*roe,x_real_store/wlh,kvec*roe,x_imag_store/wlh);
h(1).XLabel.String = 'k\rho_e';
h(1).YLabel.String = '\omega/\omega_{LH}';
h(2).YLabel.String = '\gamma/\omega_{LH}';
s1 = sprintf('B = %g nT \nn = %g cc \nTe = %g eV \nTi = %g eV \nbeta = %g \nbeta_i = %g \nbeta_e = %g\nvti = %g km/s\nroe = %.1f km',B*1e9,n*1e-6,Te,Ti,beta,betai,betae,vti*1e-3,roe*1e-3);
%s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \n\n',Ln*1e-3,LB*1e-3,LT*1e-3);
%s3 = sprintf('vn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \n',vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3);
s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \nvn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \nE = %.0f mV/m \n',Ln*1e-3,LB*1e-3,LT*1e-3,vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3,E*1e3);
s3 = sprintf('k_{max}*roe = %g \nwi_{max}/wlh = %.3f \nwr_{max}/wlh = %.3f \nvph_{max} = %.0f km/s',kmax*roe,wimax/wlh,wrmax/wlh,vph*1e-3);

text(h(1).XLim(1),h(1).YLim(2),s1,'verticalalignment','top','horizontalalignment','left')
text(mean(h(1).XLim),h(1).YLim(2),s2,'verticalalignment','top','horizontalalignment','left')
text(h(1).XLim(2),h(1).YLim(2),s3,'verticalalignment','top','horizontalalignment','right')

h(2).YGrid = 'on';
h(2).YMinorGrid = 'on';


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
end