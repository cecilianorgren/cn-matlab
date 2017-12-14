if 0
%% 2D matrix with Te and E
n = 0.06e6;
ne = n;
ni = n;
Te_vec = 500:100:2000; Te_K_vec = Te_vec*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti_vec = 7000; Ti_K_vec = Ti_vec*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
B = 10e-9; 

E_vec = (-50:5:-10)*1e-3;
LnLT = 0.1; % Ln/LT, Ln generally smaller then LT

nTe = numel(Te_vec);
nE = numel(E_vec);

nk = 200;
kvec = linspace(1e-1,2,nk)/roe;


wr_all = nan(nE,nTe);
wi_all = nan(nE,nTe);
k_all = nan(nE,nTe);
vph_all = nan(nE,nTe);

for iE = 1:nE
  for iTe = 1:nTe
    Te = Te_vec(iTe);
    Te_K = Te_K_vec(iTe);
    E = E_vec(iE);
    mms_2017Jul11.solver_Davidson1977
    wr_all(iE,iTe) = wrmax;
    wi_all(iE,iTe) = wimax;
    k_all(iE,iTe) = kmax;
    vph_all(iE,iTe) = vphmax;
  end
end

roe_all = sqrt(2*qe*Te_vec/me)/(qe*B/me);
roe_mat = repmat(roe_all,nE,1);

%% Plot
figure(104)
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;
if 1
  hca = h(isub); isub = isub + 1;  
  s1 = sprintf('B = %g nT \nn = %g cc \nTi = %g eV \nbeta_i = %g \nvti = %g km/s',B*1e9,n*1e-6,Ti,betai,vti*1e-3);
  %s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \nvn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \nE = %.0f mV/m \n',Ln*1e-3,LB*1e-3,LT*1e-3,vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3,E*1e3);  

  if exist('ht1','var'); delete(ht1); end
  %if exist('ht2','var'); delete(ht2); end
  ht1 = text(hca,hca.XLim(1),hca.YLim(2),s1,'verticalalignment','top','horizontalalignment','left');
  %ht2 = text(hca,hca.XLim(2),hca.YLim(2),s2,'verticalalignment','top','horizontalalignment','left');  
  hca.Visible = 'off';
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,k_all*1e3); 
  %shading(hca,'flat');
  hca.Title.String = 'k at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '1/km';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,k_all.*roe_mat); 
  %shading(hca,'flat');
  hca.Title.String = 'k*roe at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'k*roe';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wr_all/wlh); 
  %shading(hca,'flat');
  hca.Title.String = 'w_r at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wi_all/wlh); 
  %shading(hca,'flat');
  hca.Title.String = 'w_r at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wr_all./k_all*1e-3); 
  %shading(hca,'flat');
  hca.Title.String = 'v_{ph} at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'km/s';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
end
%% run for  array of Te and E, 3rd dim is k
units = irf_units;
n = 0.07e6;
ne = n;
ni = n;
Te_vec = 200:200:3000; Te_K_vec = Te_vec*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti = 7000; Ti_K = Ti*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
B = 10e-9; 

E_vec = (-50:10:-10)*1e-3;
E_vec = (-50)*1e-3;
LnLT = 0.1; % Ln/LT, Ln generally smaller then LT

nTe = numel(Te_vec);
nE = numel(E_vec);

nk = 250;
k_min= 1e-1;
k_max = 5;

wr_allmax = nan(nE,nTe);
wi_allmax = nan(nE,nTe);
k_allmax = nan(nE,nTe);
vph_allmax = nan(nE,nTe);

wr_all = nan(nE,nTe,nk);
wi_all = nan(nE,nTe,nk);
k_all = nan(nE,nTe,nk);
vph_all = nan(nE,nTe,nk);

for iE = 1:nE
  for iTe = 1:nTe
    Te = Te_vec(iTe);
    Te_K = Te_K_vec(iTe);
    E = E_vec(iE);
    guess_previous_k = 1;
    xguess = 0;
    mms_2017Jul11.solver_Davidson1977
    if iTe > 1
      % Compare subsequent solutions to see if they differ too much. If
      % they do, use the previous solution as initial guess (for each k)
      % for the new solution.
      wr_diff = tocolumn(squeeze(wr_all(iE,iTe-1,:))) - tocolumn(wr_store);
      wi_diff = tocolumn(squeeze(wi_all(iE,iTe-1,:))) - tocolumn(wi_store);
      x_diff = mean(abs(wi_diff));
      figure(99); plotyy(kvec,[tocolumn(squeeze(wr_all(iE,iTe-1,:))) tocolumn(squeeze(wr_store))],...
                         kvec,[tocolumn(squeeze(wi_all(iE,iTe-1,:))) tocolumn(squeeze(wi_store))]);
               title(gca,sprintf('xdiff = %g, Te = %g & %g',x_diff,Te,Te_vec(iTe-1)))
               drawnow
      while x_diff > 0.5*wlh
        fprintf('Possibly bad solution wr_diff = %g. Rerun with previous solution as initial guess.\n',x_diff)
        guess_previous_k = 0; % use entire previous solution instead      
        mms_2017Jul11.solver_Davidson1977
        wr_diff = tocolumn(squeeze(wr_all(iE,iTe-1,:))) - tocolumn(wr_store);
        wi_diff = tocolumn(squeeze(wi_all(iE,iTe-1,:))) - tocolumn(wi_store);
      figure(99); plotyy(kvec,[tocolumn(squeeze(wr_all(iE,iTe-1,:))) tocolumn(squeeze(wr_store))],...
                         kvec,[tocolumn(squeeze(wi_all(iE,iTe-1,:))) tocolumn(squeeze(wi_store))])
                       title(gca,sprintf('xdiff = %g, Te = %g & %g',x_diff,Te,Te_vec(iTe-1)))
                       drawnow
        %abs(mean(wr_diff))
        x_diff = mean(abs(wi_diff));
      end
    end
    if 0%iE > 1
      % Compare subsequent solutions to see if they differ too much. If
      % they do, use the previous solution as initial guess (for each k)
      % for the new solution.
      wr_diff = tocolumn(squeeze(wr_all(iE-1,iTe,:))) - tocolumn(wr_store);
      wi_diff = tocolumn(squeeze(wi_all(iE-1,iTe,:))) - tocolumn(wi_store);
      x_diff = mean(abs(wi_diff));
      figure(99); plotyy(kvec,[tocolumn(squeeze(wr_all(iE-1,iTe,:))) tocolumn(squeeze(wr_store))],...
                         kvec,[tocolumn(squeeze(wi_all(iE-1,iTe,:))) tocolumn(squeeze(wi_store))]);
               title(gca,sprintf('xdiff = %g, E = %g & %g',x_diff,E,E_vec(iE-1)))
               drawnow
      while x_diff > 0.5*wlh
        fprintf('Possibly bad solution wr_diff = %g. Rerun with previous solution as initial guess.\n',x_diff)
        guess_previous_k = 0; % use entire previous solution instead      
        mms_2017Jul11.solver_Davidson1977
        wr_diff = tocolumn(squeeze(wr_all(iE-1,iTe,:))) - tocolumn(wr_store);
        wi_diff = tocolumn(squeeze(wi_all(iE-1,iTe,:))) - tocolumn(wi_store);
      figure(99); plotyy(kvec,[tocolumn(squeeze(wr_all(iE-1,iTe,:))) tocolumn(squeeze(wr_store))],...
                         kvec,[tocolumn(squeeze(wi_all(iE-1,iTe,:))) tocolumn(squeeze(wi_store))])
                       title(gca,sprintf('xdiff = %g, E = %g & %g',x_diff,E,E_vec(iE-1)))
                       drawnow
        %abs(mean(wr_diff))
        x_diff = mean(abs(wi_diff));
      end
    end
    wr_allmax(iE,iTe) = wrmax;
    wi_allmax(iE,iTe) = wimax;
    k_allmax(iE,iTe) = kmax;
    vph_allmax(iE,iTe) = vphmax;
    
    wr_all(iE,iTe,:) = wr_store;
    wi_all(iE,iTe,:) = wi_store;
    k_all(iE,iTe,:) = kvec;
    vph_all(iE,iTe,:) = wr_store./kvec;
  end
end

roe_all = sqrt(2*qe*Te_vec/me)/(qe*B/me);
roe_mat = repmat(roe_all,nE,nk);
disp('ready!')

%% Plot range of Te
iE = 1;
figure(104)
nrows = 4;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;
if 1 % em coupling parameter 'd_gradB' vs k*roe
  hca = h(isub); isub = isub + 1;
  iE = 1; 
  iTe = 1:nTe;
  Pe_vec = Te_K_vec*units.kB*n; 
  betae_vec = Pe_vec/PB;
  betae_all = repmat(tocolumn(betae_vec),1,nk);
  kroe_all = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);  
  pcolor(hca,kroe_all,Te_vec(iTe),betae_all./kroe_all);   
  hold(hca,'on')
  gridx = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  shading(hca,'flat');
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = '\delta_{em} = \beta_e/k^2\rho_e^2';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k\rho_e';    
  hca.CLim(1) = 0;
end
if 1 % em coupling parameter 'd_em' vs k*roe
  hca = h(isub); isub = isub + 1;
  iE = 1; 
  iTe = 1:nTe;
  Pe_vec = Te_K_vec*units.kB*n; 
  betae_vec = Pe_vec/PB;
  betae_all = repmat(tocolumn(betae_vec),1,nk);
  Te_all = repmat(tocolumn(Te_vec),1,nk);
  kroe_all = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);  
  pcolor(hca,kroe_all,Te_vec(iTe),-0.5*betae_all.*(1+Te_all/Ti));   
  hold(hca,'on')
  gridx = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hca.Title.String = '\delta_{\nabla B} = v_B/vE = -0.5\beta_e(1+T_e/T_i)';
  hcb.YLabel.String = '\delta_{\nabla B} = -0.5\beta_e(1+T_e/T_i)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe(T_e)';    
  hca.CLim = [min(min(-0.5*betae_all.*(1+Te_all/Ti))) 0];
end
if 1 % input info
  hca = h(isub); isub = isub + 1;  
  s1 = sprintf('B = %g nT \nn = %g cc \nTi = %g eV \nbeta_i = %g \nvti = %g km/s \nE = %g mV/m \nV_E = %.0f km/s \nV_E/vti = %.2f \nwpe/wce = %.2f \n(wpe/wce)^2 = %.2f',B*1e9,n*1e-6,Ti,betai,vti*1e-3,E_vec(iE)*1e3,E_vec(iE)/B*1e-3,E_vec(iE)/B/vti,wpe/wce,(wpe/wce)^2);
  %s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \nvn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \nE = %.0f mV/m \n',Ln*1e-3,LB*1e-3,LT*1e-3,vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3,E*1e3);  

  if exist('ht1','var'); delete(ht1); end
  %if exist('ht2','var'); delete(ht2); end
  ht1 = text(hca,hca.XLim(1),hca.YLim(2),s1,'verticalalignment','top','horizontalalignment','left');
  %ht2 = text(hca,hca.XLim(2),hca.YLim(2),s2,'verticalalignment','top','horizontalalignment','left');  
  hca.Visible = 'off';
end
if 0 % empty
  hca = h(isub); isub = isub + 1;  
  hca.Visible = 'off';
end
if 0 % empty
  hca = h(isub); isub = isub + 1;  
  hca.Visible = 'off';
end
if 1 % wi vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh);   
  hold(hca,'on')
  gridx = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe(T_e)';    
end
if 1 % wi vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(k_all(iE,:,:))*1e3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % wi vs lambda = 2pi/k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(2*pi/k_all(iE,:,:))*1e-3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end
if 1 % wr vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');  
  hold(hca,'on')
  gridx = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe(T_e)';    
end
if 1 % wr vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(k_all(iE,:,:))*1e3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % wr vs lambda = 2pi/k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(2*pi/k_all(iE,:,:))*1e-3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end
if 1 % vph vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hold(hca,'on')
  gridx = squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk);
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe';    
end
if 1 % vph vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(k_all(iE,:,:))*1e3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % vph vs lambda
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hold(hca,'on')  
  gridx = squeeze(2*pi/k_all(iE,:,:))*1e-3;
  gridy = repmat(tocolumn(Te_vec),1,nk);  
  contour(hca,gridx,gridy,squeeze(wi_all(iE,iTe,:))/wlh,'k'); 
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end

cmap = cn.cmap('white_parula');

for ip = 4:npanels
  colormap(h(ip),cmap)
  h(ip).CLim = 0.4*[-1 1]*max(abs(h(ip).CLim));
end

%% Plot range of E
iTe = 1;
figure(104)
nrows = 4;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;
if 1 % input info
  hca = h(isub); isub = isub + 1;  
  s1 = sprintf('B = %g nT \nn = %g cc \nTi = %g eV \nbeta_i = %g \nvti = %g km/s \nE = %g mV/m \nV_E = %.0f km/s \nV_E/vti = %.2f \nwpe/wce = %.2f \n(wpe/wce)^2 = %.2f',B*1e9,n*1e-6,Ti,betai,vti*1e-3,E_vec(iE)*1e3,E_vec(iE)/B*1e-3,E_vec(iE)/B/vti,wpe/wce,(wpe/wce)^2);
  %s2 = sprintf('Ln = %g km \nLB = %g km \nLT = %g km \nvn = %.0f km/s \nvT = %.0f km/s \nvB = %.0f km/s \nvE = %.0f km/s \nE = %.0f mV/m \n',Ln*1e-3,LB*1e-3,LT*1e-3,vn*1e-3,vT*1e-3,vB*1e-3,vE*1e-3,E*1e3);  

  if exist('ht1','var'); delete(ht1); end
  %if exist('ht2','var'); delete(ht2); end
  ht1 = text(hca,hca.XLim(1),hca.YLim(2),s1,'verticalalignment','top','horizontalalignment','left');
  %ht2 = text(hca,hca.XLim(2),hca.YLim(2),s2,'verticalalignment','top','horizontalalignment','left');  
  hca.Visible = 'off';
end
if 1 % empty
  hca = h(isub); isub = isub + 1;  
  hca.Visible = 'off';
end
if 1 % empty
  hca = h(isub); isub = isub + 1;  
  hca.Visible = 'off';
end
if 1 % wi vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe(T_e)';    
end
if 1 % wi vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % wi vs lambda = 2pi/k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(wi_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end
if 1 % wr vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe(T_e)';    
end
if 1 % wr vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % wr vs lambda = 2pi/k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(wr_all(iE,iTe,:))/wlh); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end
if 1 % vph vs kroe
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:)).*repmat(tocolumn(roe_all),1,nk),Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k*roe';    
end
if 1 % vph vs k
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,squeeze(k_all(iE,iTe,:))*1e3,Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = 'k (1/km)';    
end
if 1 % vph vs lambda
  iE = 1; 
  iTe = 1:nTe;
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,2*pi./squeeze(k_all(iE,iTe,:))*1e-3,Te_vec(iTe),squeeze(vph_all(iE,iTe,:))*1e-3); 
  shading(hca,'flat');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ph} (km/s)';
  hca.YLabel.String = 'T_e (eV)';
  hca.XLabel.String = '\lambda (km)';    
end

cmap = cn.cmap('white_parula');

for ip = 4:npanels
  colormap(h(ip),cmap)
  h(ip).CLim = [-1 1]*max(abs(h(ip).CLim));
end
%%
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,k_all.*roe_mat); 
  %shading(hca,'flat');
  hca.Title.String = 'k*roe at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'k*roe';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wr_all/wlh); 
  %shading(hca,'flat');
  hca.Title.String = 'w_r at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_r/w_{LH}';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wi_all/wlh); 
  %shading(hca,'flat');
  hca.Title.String = 'w_r at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'w_i/w_{LH}';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,Te_vec,E_vec*1e3,wr_all./k_all*1e-3); 
  %shading(hca,'flat');
  hca.Title.String = 'v_{ph} at maximum w_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'km/s';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'T_e (eV)';    
end