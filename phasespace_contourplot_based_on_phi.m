localuser = datastore('local','user');
% Load MMS data
if 0 % load data
  %%
  ic = 1:4;
  tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
  mms.db_init('local_file_db','/Volumes/Nexus/data');
  db_info = datastore('mms_db');   
  c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
  c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
  c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
  c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
  c_eval('intEdt? = irf_integrate(gseE?par);');    
  mms.load_data_edi;  
  c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
  % Make reduced distribution
  tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.617Z');
  tintZoom = tint_phi + [-2 2];   
  strTintZoom = [irf_time(tintZoom(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tintZoom(2),'epochtt>utc_HHMMSS')];
  eint = [000 40000];
  vint = [-Inf Inf];
  eDist = ePDist1.tlim(tintZoom).elim(eint);   
  scpot = scPot1.resample(eDist);
  scpot_margin = 1.0; % keep in mind that this also affects the velocity at lower energies
  lowerelim = scpot*scpot_margin;
  eLine = dmpaB1.resample(eDist).norm;
  tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B  
end

% Plotting options
doT = 1; % otherwise plot x;

% Plasma properties
units = irf_units;
n = [0.045]*1e6;
T = [400]; T_K = T*units.eV/units.kB; % use parallel temperature
vd = [-0000]*1e3; % m/s

% EDI energy and corresponding velocity
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV
dv_edi = sqrt(2*units.e*dE_edi./units.me); % m/s

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
dv_edi_plus = sqrt(2*units.e*E_edi_plus./units.me)-v_edi; % m/s
dv_edi_minus = -sqrt(2*units.e*E_edi_minus./units.me)+v_edi; % m/s
if 0 % illustrate energy and velocity range
  
end

% Physical parameters
vt = sqrt(2*units.e*T./units.me); % m/s
wp = sqrt(n.*units.e^2./(units.me*units.eps0)); % Hz
Le = units.c./wp;
Ld = vt./wp/sqrt(2);

% Wave properties
% potential
% x can be vph * t if w e want to comapre to time series
% spatial grid (can also be time series, i.e. use directly)

% eh model
% observed/measured properties (konrad)
data_tmp = load('/Users/cno062/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat');
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(velocity);
c_eval('obs_t0_epoch_mms? = EH_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')

%
% Potential from observed E
tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.620Z');
%tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.700Z');
vph = -9000e3; % m/s, representative phase velocity
ffilt = 20; % Hz
phi_shift = 150; % to keep potential > 0
c_eval('Etoint? = gseE?par.filt(ffilt,0,[],3);');    
c_eval('intEdt? = irf_integrate(Etoint?.tlim(tint_phi));');    
c_eval('phi? = intEdt?*vph*1e-3;')
minpeakdistance = 150;
c_eval('[PKS?,LOCS?,W?] = findpeaks(-phi?.data,''MinPeakDistance'',minpeakdistance);')

c_eval('phi?_detrend = phi?; phi?_detrend.data = detrend(phi?_detrend.data,''linear'',LOCS?);')
c_eval('phi?_detrend_shift = phi?_detrend + phi_shift;')
c_eval('phi?_detrend_shift.data(phi?_detrend_shift.data<0) = 0;')

if 0 % plot different stages of phi: original from filtered E, detrend, shift
figure(31)
h = irf_plot(7); isub = 1;
irf_plot(h(isub),{phi1,phi2,phi3,phi4},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi1_detrend,phi2_detrend,phi3_detrend,phi4_detrend},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi1_detrend_shift,phi2_detrend_shift,phi3_detrend_shift,phi4_detrend_shift},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi1,phi1_detrend,phi1_detrend_shift,phi1-phi1_detrend},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi2,phi2_detrend,phi2_detrend_shift,phi2-phi2_detrend},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi3,phi3_detrend,phi3_detrend_shift,phi3-phi3_detrend},'comp'); isub = isub + 1;
irf_plot(h(isub),{phi4,phi4_detrend,phi4_detrend_shift,phi4-phi4_detrend},'comp'); isub = isub + 1;
irf_zoom(h,'x',tint_phi)
end

t0 = tint_phi(1);
x_vec = phi1.time - t0; % seconds
nx = numel(x_vec);

mms_id = 1; % chose potentiual between spacecraft
c_eval('phi_vec = phi?_detrend_shift.data;',mms_id)
c_eval('phi_vec_obs = obs_phi?.data;',mms_id)
c_eval('epar_vec = Etoint?.tlim(tint_phi).data;',mms_id)
c_eval('vph_vec_obs = obs_vph?.data;',mms_id)
c_eval('x_vec_obs = obs_phi?.time-t0;',mms_id)
c_eval('LOCS = LOCS?;',mms_id)
c_eval('PKS = PKS?;',mms_id)

option_vph = 'interp linear piecewise';
%option_vph = 'constant';
switch option_vph
  case 'gradual'
    vph_vec = vph*(1+0.5*x_vec/max(x_vec));
  case 'constant'
    vph_vec = repmat(vph,1,nx);
  case 'interp linear piecewise' % piecewise linear interpolation of vph
    c_eval('tmp_time = obs_t0_epoch_mms?-t0;',mms_id)
    tmp_time = torow([tint_phi(1)-t0; tmp_time]);
    c_eval('tmp_data = obs_vph?.data*1e3;',mms_id) % m/s
    tmp_data = torow([vph; tmp_data]);        
    vph_vec = interp_linear_piecewise(tmp_data,tmp_time,x_vec);
    vph_vec = smooth(vph_vec,numel(x_vec)/50);
    %plot(x_vec,vph_vec,'.',tmp_time,tmp_data,'*')
end
% adjust phi incase vph is variable (vph is used to get phi: phi = eint*vph)
phi_vec = phi_vec.*reshape(vph_vec,size(phi_vec))/vph;

% velocity grid
vmax = 2*vt;
nv = 1500;
v = linspace(-vmax,vmax,nv);
dv = v(2) - v(1);
t_to_x = -vph;
if doT  
  dx_vec = x_vec(2) - x_vec(1);
  dx = dx_vec*t_to_x;
else
  dx_vec = x_vec(2) - x_vec(1);
  dx = x_vec(2) - x_vec(1);
end

[X,V] = meshgrid(x_vec,v);
VPH = repmat(torow(vph_vec),nv,1);
PHI = repmat(torow(phi_vec),nv,1);
 
beta = -5; % decides phase space depletion inside EH
F = maxwellian_phi(V,n,vt,vd,PHI,VPH,beta); % distribution function, from Schamel 1986
Fdv = F*dv;
FV = F.*V; % (s1/m4)*(m1/s1) = (1/m3)
FV2 = F.*V.^2; % (s1/m4)*(m1/s1)^2 = (1/m2s)
FVdv = FV*dv; % (s1/m4)*(m1/s1)^2 = (1/m2s)
FV2dv = FV2*dv; % (s1/m4)*(m1/s1)^3 = (1/ms2)

% 'integrals'
sumFdv = nansum(Fdv,1); % (s1/m4)*(m1/s1) = (1/m3) 
sumFVdv = nansum(FVdv,1); % (s1/m4)*(m1/s1)*(m1/s1) = (1/m2s)
sumFV2dv = nansum(FV2dv,1); % (s1/m4)*(m1/s1)*(m1/s1) = (1/m2s)

%
% charge density from observed phi
x_vec_diff1 = x_vec(1:end-1)+0.5*dx_vec;
x_vec_diff2 = x_vec(2:end-1);
efield_from_obs_phi = -diff(phi_vec,1)/dx;
density_diff_from_obs_phi = diff(phi_vec,2)*units.eps0/units.e/dx/dx;
electron_density_from_obs_phi = n*units.e-charge_density_from_obs_phi;



% phi = phi0*exp(-x^2/2/l^2), 
% dphi/dx = (-2*x/2/l^2)*phi0*exp(-x^2/2/l^2)
% d2phi/dx2 = (-1/l^2)*phi0*exp(-x^2/2/l^2) + (-2*x/2/l^2)^2*phi0*exp(-x^2/2/l^2) 
%           = [(-1/l^2) + (-x/l^2)^2]*phi0*exp(-x^2/2/l^2) 
%           = -(1/l^2)[1 + (x/l)^2]*phi0*exp(-x^2/2/l^2) 
%           == -(e/eps0)*(ni-ne)
% @ x = 0: 
% d2phi/dx2 = -(1/l^2)*phi0 = -(e/eps0)*(ni-ne)
%   => (ni-ne) = eps0/e*(1/l^2)*phi0
% potential from charge density
dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc

% Poissons equation, solve to find proper beta
% integrate f: sumFdv
model_density = sumFdv; % (s1/m4)*(m1/s1) = (1/m3) 
model_charge_density = units.e*(n-model_density); % rho = e(ni-ne)
model_charge_density = units.e*(mean(model_density)-model_density); % rho = e(ni-ne)
model_efield = cumsum(  model_charge_density-0*mean(model_charge_density))/units.eps0*dx; % V/m
model_efield = cumtrapz(model_charge_density-0*mean(model_charge_density))/units.eps0*dx; % V/m
model_phi = -cumsum(model_efield)*dx;

%efield = sign(vph)*(1/units.eps0)*tocolumn(cumtrapz(x_vec,charge_density_prep))*t_to_x;
%potential = -1*tocolumn(cumtrapz(x_vec,efield))*t_to_x;
%potential = detrend(potential,'linear');
%potential_fit = polyfit(x_vec,potential,3);
%potential_dc = polyval(potential_fit,x_vec);
if 1 % plot derivation of beta
  figure(32)
  clear h
  nrows = 7;
  ncols = 1;
  npanels = nrows*ncols;
  for ip = 1:npanels
    h(ip) = subplot(nrows,ncols,ip);
  end
  isub = 1;
  if 1 % efield
    hca = h(isub); isub = isub + 1; 
    plot(hca,x_vec,epar_vec,efield_from_obs_phi_x_vec,efield_from_obs_phi*1e3)
    irf_legend(hca,{'E_{obs}'},[0.01 0.99])
    hca.YLabel.String = {'E','mV/m'};
  end
  if 1 % density 
    hca = h(isub); isub = isub + 1; 
    plotyy(hca,x_vec,model_density*1e-6,x_vec_diff2,density_diff_from_obs_phi*1e-6+mean(model_density)*1e-6)
    irf_legend(hca,{'n_{e,mod}','n_e-n_i'},[0.01 0.99])
    hca.YLabel.String = {'n','cm^{-3}'};
  end
  hca = h(isub); isub = isub + 1; 
  plotyy(hca,x_vec,density,charge_density_from_obs_phi_x_vec,charge_density_from_obs_phi*1e3)
  hca = h(isub); isub = isub + 1; 
  plotyy(hca,x_vec,density,charge_density_from_obs_phi_x_vec,electron_density_from_obs_phi*1e3)
  hca = h(isub); isub = isub + 1; 
  plotyy(hca,x_vec,model_density,x_vec,model_charge_density)  
  hca = h(isub); isub = isub + 1; 
  plotyy(hca,x_vec,epar_vec,x_vec,detrend(model_efield*1e3,'linear',1:40:nx))  
  hca = h(isub); isub = isub + 1; 
  plotyy(hca,x_vec,phi_vec,x_vec,detrend(model_phi))  
end
%%
% flux at edi energy/velocity interval
vind_edi_0 = intersect(find(v>v_edi-dv_edi),find(v<v_edi+dv_edi));
vind_edi_180 = intersect(find(v>-v_edi-dv_edi),find(v<-v_edi+dv_edi));
FVdv_edi_0 = nansum(FVdv(vind_edi_0,:));
FVdv_edi_180 = nansum(FVdv(vind_edi_180,:));


% plot
figure(33)
clear h;
nrows = 7;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

c_eval('edi_tshift = 0*0.5*(flux180_mms?.time(2)-flux180_mms?.time(1));',ic)
if 0 % PHI, pcolor
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_vec,V*1e-6,PHI);
  shading(hca,'flat')
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '\phi';
end
if 1 % Epar, plot
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,epar_vec);  
  hca.YLabel.String = 'E_{||} (mV/m)';  
end
if 0 % Epar, model
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,efield*1e3);  
  irf_legend(hca,{'from model'},[0.01 0.99],'color',[0 0 0]);
  hca.YLabel.String = 'E (mV/m)';  
end
if 1 % PHI, plot
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,phi_vec);  
  if 1
    hold(hca,'on')
    plot(hca,x_vec(LOCS),phi_vec(LOCS),'*')
    hold(hca,'off')
  end
  irf_legend(hca,{sprintf('v_{ph}=%gx10^3 km/s, phi_{shift}=%g V, f_{filt,E}=%g Hz',vph*1e-6,phi_shift,ffilt)},[0.01 0.99],'color',[0 0 0]);
  hca.YLabel.String = '\phi (V)';  
end
if 0 % PHI, model
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,potential);  
  irf_legend(hca,{'from model'},[0.01 0.99],'color',[0 0 0]);
  hca.YLabel.String = '\phi (V)';  
end
if 1 % F
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,(F));
  shading(hca,'flat')
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f (s^1/m^4)';
  colormap(hca,cn.cmap('white_blue'))
  %colormap(hca,cn.cmap('blue_white'))
  
  if 1 % EDI energies
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; [0.9290    0.6940    0.1250];
    hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end
    hlines = plot(hca,...
      x_vec([1 end]),(v_edi-dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(v_edi+dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi-dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi+dv_edi)*1e-6*[1 1],...
      'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
    irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 1 % model phase velocity
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
    hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
    irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 1 % observed phase velocity
    hold(hca,'on')
    hlines = plot(hca,x_vec_obs,vph_vec_obs*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
    irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end  
  
  irf_legend(hca,{sprintf('T_{bg}=%g eV, n_{bg}=%g cc, beta_{Schamel}=%g',T,n*1e-6,beta)},[0.99 0.99],'color',hlines(1).Color);    
end
if 0 % e psd vpar
  hca = h(isub); isub = isub + 1;
  specrec = ef1D.tlim(tint_phi).specrec('velocity_1D','10^3 km/s');
  time_fred = ef1D.tlim(tint_phi).time-tint_phi(1);
  imagesc(hca,time_fred,specrec.f(1,:),specrec.p');
  hcb = colorbar('peer',hca);
  colormap(cn.cmap('white_blue'))
  hca.CLim = [0 3]*1e-3; 
  
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % log10(F)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,log10(F));
  shading(hca,'flat')
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}f (s^1/m^4)';
  colormap(hca,cn.cmap('white_blue'))
  %colormap(hca,cn.cmap('blue_white'))
  hold(hca,'on')
  line_color = [0.5 0.5 0.5]; [0.9290    0.6940    0.1250];
  %hvph = plot(hca,x_vec,vph_vec*1e-6,'Color',line_color,'LineWidth',1.5);
  %irf_legend(hca,{'v_{ph}'},[0.01 0.1],'color',line_color);
  hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
  hlines(1).LineStyle = '--'; hlines(1).Color = [0 0 0];
  hlines(2).LineStyle = '--'; hlines(2).Color = [0 0 0];
  irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);  
  hold(hca,'off')
end
if 0 % e psd vpar log 10
  hca = h(isub); isub = isub + 1;
  specrec = ef1D.tlim(tint_phi).specrec('velocity_1D','10^3 km/s');
  time_fred = ef1D.tlim(tint_phi).time-tint_phi(1);
  imagesc(hca,time_fred,specrec.f(1,:),log10(specrec.p)');
  hcb = colorbar('peer',hca);
  colormap(cn.cmap('white_blue'))
  hca.CLim = [-5.5 -3.5];   
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % sum F (density)
  hca = h(isub); isub = isub + 1;  
  plot(hca,x_vec,sumFdv*1e-6); % 1e-6 from 1/m3 > 1/cm3
  hold(hca,'on')
  h_nmean = plot(hca,x_vec([1 end]),mean(sumFdv*1e-6)*[1 1]); % 1e-6 from 1/m3 > 1/cm3
  irf_legend(hca,{'<n>'},[0.01 0.99],'color',h_nmean(1).Color);  
  hold(hca,'off')
  hca.YLabel.String = 'n (cm^{-3})';
  %yticks = hca.YTick;
  %yticklabels = num2str(tocolumn(yticks),'%.3f');
  %hca.YTickLabels = yticklabels;
end
if 1 % flux: F*v
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,FV);
  shading(hca,'flat')
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f*v (1/m^3)'; 
  hca.CLim = hca.CLim(2)*[-1 1]; 
  hcb.YLim = hca.CLim;
  if 1 % EDI energies
    hold(hca,'on')
    hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
    hlines = plot(hca,...
      x_vec([1 end]),(v_edi-dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(v_edi+dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi-dv_edi)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi+dv_edi)*1e-6*[1 1],...
      'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
    irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
    hold(hca,'off')
  end
  if 1 % model phase velocity
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
    hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
    irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 1 % observed phase velocity
    hold(hca,'on')
    hlines = plot(hca,x_vec_obs,vph_vec_obs*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
    irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end  
  colormap(hca,cn.cmap('blue_red'))
end
if 0 % flux at edi energy
  hca = h(isub); isub = isub + 1;
  units_scale = 1e-4;
  plot(hca,x_vec,abs(FVdv_edi_0)*units_scale,x_vec,abs(FVdv_edi_180)*units_scale,x_vec,(abs(FVdv_edi_180)-abs(FVdv_edi_0))*units_scale);  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'0^o','180^o','180^o-0^o','location','eastoutside')
end
if 0 % flux at edi energy
  hca = h(isub); isub = isub + 1;
  units_scale = 1e-4;
  plot(hca,x_vec,abs(FVdv_edi_0)*units_scale,x_vec,abs(FVdv_edi_180)*units_scale);  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'0^o','180^o','location','eastoutside')  
  text(hca,0,hca.YLim(2),'model','verticalalignment','top')
end
if 0 % flux measured by EDI
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-flux0.time(1);
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-flux180.time(1);
  flux180_data = flux180.data;
  nodes = 1:2;
  plot(hca,flux0_time,mean(flux0_data(:,nodes),2),flux180_time,mean(flux180_data(:,nodes),2));  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'0^o','180^o','location','eastoutside')
  text(hca,0,hca.YLim(2),'EDI','verticalalignment','top')
end
if 0 % cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-flux0.time(1);
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-flux180.time(1);
  flux180_data = flux180.data;
  nodes = 1:2;
  units_scale = 1e-4; % m^-2 > cm^-2 
  plot(hca,x_vec,abs(FVdv_edi_180)*units_scale,flux180_time + edi_tshift,mean(flux180_data(:,nodes),2));  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'model','EDI','location','eastoutside')
  text(hca,0,hca.YLim(2),'180^o','verticalalignment','top')
end
if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-t0;
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-t0;
  flux180_data = flux180.data;
  nodes = 1:2;
  units_scale = 1e-4; % m^-2 > cm^-2 
  units_scale_2 = 1e6;
  plot(hca,x_vec,abs(FVdv_edi_180)*units_scale/units_scale_2,flux180_time + edi_tshift,mean(flux180_data(:,nodes),2)/units_scale_2);  
  %hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
  legend(hca,'model','EDI','location','eastoutside')
  text(hca,0,hca.YLim(2),'180^o','verticalalignment','top')
end
if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-t0;
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-t0;
  flux180_data = flux180.data;
  nodes = 1:2;
  units_scale = 1e-4; % m^-2 > cm^-2 
  units_scale_2 = 1e6;
  %hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  plot(hca,x_vec,abs(FVdv_edi_0)*units_scale/units_scale_2,flux0_time + edi_tshift,mean(flux0_data(:,nodes),2)/units_scale_2);  
  hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
  legend(hca,'model','EDI','location','eastoutside')
  text(hca,0,hca.YLim(2),'0^o','verticalalignment','top')
end
if 0 % SI UNITS, comparing model flux with flux measured by EDI, at 180
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-flux0.time(1);
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-flux180.time(1);
  flux180_data = flux180.data;
  nodes = 1:2;
  units_scale = 1e-4; % m^-2 > cm^-2 
  units_scale_2 = 1e9;
  plot(hca,x_vec,abs(FVdv_edi_180)/units_scale_2,flux180_time+0.5*diff(flux180_time(1:2)),mean(flux180_data(:,nodes),2)/units_scale/units_scale_2);  
  hca.YLabel.String = sprintf('flux (10^%g m^{-2}s^{-1})',log10(units_scale_2));
  legend(hca,'model','EDI','location','eastoutside')
  text(hca,0,hca.YLim(2),'180^o','verticalalignment','top')
end
if 0 % flux at edi energy
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,abs(flux_edi_0),x_vec,abs(flux_edi_180),x_vec,abs(flux_edi_180)-abs(flux_edi_0));  
  hca.YLabel.String = 'flux (...)';
  legend(hca,'0^o','180^o','180^o-0^o','location','eastoutside')
end
if 0 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
  colormap(hca,cn.cmap('blue_red'))
end
if 0 % sum FV
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,sumFV);  
  hca.YLabel.String = 'v (#)';
end

%colormap(cn.cmap('blue_red'))

axes_width = h(1).Position(3);
for ipanel = 1:npanels
  if hca.Position(3) < axes_width
    axes_width = hca.Position(3);
  end
end

for ipanel = 1:npanels
  hca = h(ipanel);
  hca.Position(3) = axes_width;
  hca.Position(4) = hca.Position(4)*1.3;
  hca.Position(3);
  hca.XLim = x_vec([1 end]);
end
for ipanel = 1:(npanels-1)
  hca = h(ipanel);
  hca.XTickLabel = [];  
end
h(end).XLabel.String = 'time (s)';
linkaxes(h,'x')
irf_plot_axis_align
h(1).YLim = [-49 49];