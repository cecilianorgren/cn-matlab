% see also
% open lic.phasespace
% open lic.phasespace_1

% Load MMS data
if 1 % load data
  %%
  ic = 1:4;
  tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
  mms.db_init('local_file_db','/Volumes/Nexus/data');
  db_info = datastore('mms_db');   
  c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
  c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
  c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
  c_eval('intEdt? = irf_integrate(gseE?par);');    
  mms.load_data_edi;  
end

% Plotting options
doT = 1; % otherwise plot x;

% Plasma properties
units = irf_units;
n = [0.01]*1e6;
T = [500]; T_K = T*units.eV/units.kB; % use parallel temperature
vd = [-2000]*1e3; % m/s

% EDI energy and corresponding velocity
E_edi = 500;
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s

% Physical parameters
vt = sqrt(2*units.e*T./units.me); % m/s
wp = sqrt(n.*units.e^2./(units.me*units.eps0)); % Hz
Le = units.c./wp;
Ld = vt./wp/sqrt(2);

% Wave properties
% potential
% x can be vph * t if w e want to comapre to time series
% spatial grid (can also be time series, i.e. use directly)
wave_potential = 'konrad';
switch wave_potential
  case 'cos'    
    phi0 = 200; % V
    wavelength = 8e3; % m
    k = 2*pi/wavelength;
    phi = @(k,x) phi0*0.5*(1+cos(k*x+pi)); % same units as phi0
    Ex = @(k,x) k*phi0*(sin(k*x+pi)); % same units as phi0
    t = 0:(1/8000/10):0.002;
    x = t*vph;
    nx = numel(x);
  case 'konrad'        
    EH_properties = load('/Users/cecilia/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat');
end

% phase velocity
wave_vph = 'intE const vph';
switch wave_vph
  case 'gradual'
    vph = 10000e3; % m/s
    vph_vec = vph*(1+1*x/max(x));
  case 'steps'
    vph = 10000e3; % m/s
    vph_vec = vph*ones(1,nx);
    vph_vec(1:64) = vph_vec(1:64)*3;
    vph_vec(65:130) = vph_vec(65:130)*0.5;
  case 'konrad'    
    mms_id = 3;
    %load('/Users/cno062/GoogleDrive/Data/Events/2017-07-06_081603/')
    EH_properties = load('/Users/cecilia/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat');
    EH_properties = EH_properties.EH_properties;
    neh = numel(EH_properties.velocity);
    ts_epar = irf.ts_scalar([],[]);
    ts_vph = irf.ts_scalar([],[]);
    ts_phi = irf.ts_scalar([],[]);
    ts_vph_ = irf.ts_scalar(EH_properties.time_mms1,EH_properties.velocity);
    for ieh = 1:neh            
      eh_T = 0.04;
      c_eval('tint_tmp? = EH_properties.time_mms!(?) + eh_T*[-0.05 0.05];',ieh);      
      c_eval('tint_tmp = tint_tmp?;',ieh)
      c_eval('tidx_tmp = gseE?par.time.tlim(tint_tmp);',mms_id);      
      c_eval('[idx_tmp,epoch_tmp] = gseE?par.time.tlim(tint_tmp);',mms_id)      
      c_eval('data_tmp = gseE?par.tlim(tint_tmp).data;',mms_id);
      c_eval('eint_tmp = gseEpardt?.tlim(tint_tmp).data;',mms_id);
      c_eval('vph_tmp = EH_properties.velocity(?)*ones(1,numel(idx_tmp));',ieh);
      phi_tmp = detrend(eint_tmp,'linear')*vph_tmp(1); % all vph_tmp are the same
      %c_eval('EH_properties(?).time_mms!(?)+0.05[-0.5 0.5]',ieh,mms_id);
      
      c_eval('ts_eh? = irf.ts_scalar(epoch_tmp,data_tmp);',ieh);      
      c_eval('ts_epar = combine(ts_epar,ts_eh?);',ieh)
      
      c_eval('ts_vph? = irf.ts_scalar(epoch_tmp,vph_tmp);',ieh);
      
      c_eval('ts_vph = combine(ts_vph,ts_vph?);',ieh)
      
      c_eval('ts_phi? = irf.ts_scalar(epoch_tmp,phi_tmp);',ieh);
      c_eval('ts_phi = combine(ts_phi,ts_phi?);',ieh)
    end
    time_all = ts_Epar.time;
    c_eval('intEdt = gseEpardt?.tlim(time_all);',mms_id)    
    ts_vph_resamp = ts_vph.resample(intEdt);
    ts_Phi_nondetrend = intEdt*ts_vph_resamp;
    ts_Phi = ts_Phi_nondetrend.filt(50,0,[],3);
  %c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
  case 'konrad2' % model based on phi0    
    mms_id = 3;
    % measured properties
    EH_properties = load('/Users/cecilia/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat');
    EH_properties = EH_properties.EH_properties;
    neh = numel(EH_properties.velocity);    
    L_pp = EH_properties.L_pp;
    potential = EH_properties.calc_Pot;
    velocity = EH_properties.velocity;
    c_eval('t0_epoch = EH_properties.time_mms?;',mms_id)

    % model
    phi_mod = @(x,l,phi0) phi0*exp(-x.^2/(2*l^2));
    E_mod = @(x,l,phi0) (1/l^2)*phi0*x.*exp(-x.^2/(2*l^2));
    dt = 0.001;
    
    % collect together
    ts_epar = irf.ts_scalar([],[]);
    ts_vph = irf.ts_scalar([],[]);
    ts_phi = irf.ts_scalar([],[]);
    
    for ieh = 1:neh
      l_tmp = L_pp(ieh)/2;
      x_tmp = l_tmp*3;
      vph_tmp = velocity(ieh);
      %t_tmp = x_tmp/vph_tmp;
      phi0_tmp = potential(ieh);
      t0_epoch_tmp = t0_epoch(ieh);
      
      nx_tmp = 40;
      x_vec_tmp = linspace(-x_tmp,x_tmp,nx_tmp); % km
      vph_vec_tmp = repmat(vph_tmp,1,nx_tmp);
      t_vec_tmp = -x_vec_tmp/vph_tmp;
      t_vec_epoch_tmp = t0_epoch_tmp + t_vec_tmp;
      
      phi_tmp = phi_mod(x_vec_tmp,l_tmp,phi0_tmp);
      e_tmp = E_mod(x_vec_tmp,l_tmp,phi0_tmp);
      ts_e_tmp = irf.ts_scalar(t_vec_epoch_tmp,e_tmp);
      ts_phi_tmp = irf.ts_scalar(t_vec_epoch_tmp,phi_tmp);
      ts_vph_tmp = irf.ts_scalar(t_vec_epoch_tmp,vph_vec_tmp);
      
      %c_eval('ts_eh? = irf.ts_scalar(epoch_tmp,data_tmp);',ieh);      
      %c_eval('ts_epar = combine(ts_epar,ts_eh?);',ieh)
            
      ts_vph = combine(ts_vph,ts_vph_tmp);            
      ts_phi = combine(ts_phi,ts_phi_tmp);
      ts_epar = combine(ts_epar,ts_e_tmp);               
    end    
    
    % interp vectors    
    timeline = gseE1par.tlim(ts_phi.time).time;
    ts_vph = ts_vph.resample(timeline);
    ts_phi = ts_phi.resample(timeline);
    ts_epar = ts_epar.resample(timeline);
    
    vph_vec = ts_vph.data*1e3; % m/s
    phi_vec = ts_phi.data;
    
    if doT, x_vec = timeline - timeline(1);
    else x_vec = (timeline - timeline(1))*mean(ts_vph.data);
    end
    if 0 % plot
      c_eval('h = irf_plot({gseE?par,-1*intEdt,ts_epar,ts_phi,ts_vph});',mms_id)
      %irf_plot(h(end),ts_vph_,'*')
      irf_zoom(h,'x',ts_Epar.time);
      c_eval('hpl(!) = irf_pl_mark(h(!),tint_tmp?);',1:neh,1:numel(h))
    end
  case 'intE const vph'
    vph = -10000e3; % m/s
    tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.617Z');
    ffilt = 20;
    phi_shift = 200;
    c_eval('Etoint? = gseE?par.filt(ffilt,0,[],3);');    
    c_eval('intEdt? = irf_integrate(Etoint?.tlim(tint_phi));');    
    c_eval('phi? = intEdt?*vph*1e-3+phi_shift;')
    c_eval('phi?.data(phi?.data<0) = 0;')
    %c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
    x_vec = phi1.time - phi1.time(1);
    nx = numel(x_vec);
    vph_vec = repmat(vph,1,nx);
    phi_vec = phi1.data;
    ic = 1;
  case 'constant'
    vph_vec = vph*ones(1,nx);
end

%%


%E_vph = units.me*vph^2/2/units.eV;
%

%v_trap = sqrt(2*units.e*phi(k,x)/units.me);

% 1D Maxwellian
%f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2);
%f = @(v,n,vt,vd,phi,vph) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2 -(v-vph).^2./(2*units.e*phi/units.me).^2);
%f = @(v,n,vt,vd,phi,vph) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2 + 2*units.e*phi/units.me/vt^2);
%f = @(v,n,vt,vd,phi,vph) n*(1/pi/vt^2)^(1/2)*exp(-(sqrt(v.^2 - 2*units.e*phi/units.me) + vph).^2/vt^2);

% velocity grid
% treat this as the initial velocity of a given particle
vmax = 2*vt;
nv = 1500;
v = linspace(-vmax,vmax,nv);
Ek = units.m*v.^2/2; % kinetic energy of particle

%v_waveframe = linspace(-vmax,vmax,nv)-vph;
%Ek_waveframe = units.m*v_waveframe.^2/2; % kinetic energy of particle

[X,V] = meshgrid(x_vec,v);
%E_J = units.me*V.^2/2; % J
%E = E_J/units.e; % res=Units.me*data.^2*10^6/2/Units.eV;
VPH = repmat(torow(vph_vec),nv,1);
PHI = repmat(torow(phi_vec),nv,1);
  
%F = f(V,n,vt,vd,PHI,vph);
%F = maxwellian_phi(V-vph,n,vt,vd,PHI,vph);
beta = -5; % decides phase space depletion inside EH
F = maxwellian_phi(V-VPH+0*vd,n,vt,vd,PHI,VPH-1*vd,beta);
FV = F.*V; % (s1/m4)*(m1/s1) = (1/m3)
FV2 = F.*V.^2; % (s1/m4)*(m1/s1)^2 = (1/m2s)
dv = v(2) - v(1);
FVdv = FV*dv;
Fdv = F*dv;

% counts
sumF = nansum(F,1); % (s1/m4)
sumFdv = nansum(Fdv,1); % (s1/m4)*(m1/s1) = (1/m3) 
sumFV = nansum(FV,1); % (s1/m4)*(m1/s1) = (1/m3)
sumFVdv = nansum(FVdv,1); % (s1/m4)*(m1/s1)*(m1/s1) = (1/m2s)

% flux at edi energy/velocity
vind_edi_0 = find(abs(v-v_edi)==min(abs(v-v_edi)));
vind_edi_180 = find(abs(v+v_edi)==min(abs(v+v_edi)));
flux_edi_0 = FV2(vind_edi_0,:);
flux_edi_180 = FV2(vind_edi_180,:);
FVdv_edi_0 = FV2(vind_edi_0,:);
FVdv_edi_180 = FV2(vind_edi_180,:);

% plot
clear h;
nrows = 7;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 0 % PHI, pcolor
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_vec,V*1e-6,PHI);
  shading(hca,'flat')
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '\phi';
end
if 1 % PHI, plot
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,phi_vec);  
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
  hold(hca,'on')
  line_color = [0.5 0.5 0.5]; [0.9290    0.6940    0.1250];
  hvph = plot(hca,x_vec,vph_vec*1e-6,'Color',line_color,'LineWidth',1.5);
  irf_legend(hca,{'v_{ph}'},[0.01 0.1],'color',line_color);
  hold(hca,'off')
end
if 1 % sum F
  hca = h(isub); isub = isub + 1;  
  plot(hca,x_vec,sumFdv*1e-6); % 1e-6 from 1/m3 > 1/cm3
  hca.YLabel.String = 'n (cm^{-3})';
  yticks = hca.YTick;
  yticklabels = num2str(tocolumn(yticks),'%.3f');
  hca.YTickLabels = yticklabels;
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
  hold(hca,'on')
  hlines = plot(hca,x([1 end]),v_edi*1e-6*[1 1],x([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
  hold(hca,'off')
  colormap(hca,cn.cmap('blue_red'))
end
if 0 % flux at edi energy
  hca = h(isub); isub = isub + 1;
  units_scale = 1e-4;
  plot(hca,x_vec,abs(FVdv_edi_0)*units_scale,x_vec,abs(FVdv_edi_180)*units_scale,x_vec,(abs(FVdv_edi_180)-abs(FVdv_edi_0))*units_scale);  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'0^o','180^o','180^o-0^o','location','eastoutside')
end
if 1 % flux at edi energy
  hca = h(isub); isub = isub + 1;
  units_scale = 1e-4;
  plot(hca,x_vec,abs(FVdv_edi_0)*units_scale,x_vec,abs(FVdv_edi_180)*units_scale);  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'0^o','180^o','location','eastoutside')  
  text(hca,0,hca.YLim(2),'model','verticalalignment','top')
end
if 1 % flux measured by EDI
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
if 1 % comparing model flux with flux measured by EDI, at 180
  hca = h(isub); isub = isub + 1;
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',ic)
  flux0_time = flux0.time-flux0.time(1);
  flux0_data = flux0.data;
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',ic)
  flux180_time = flux180.time-flux180.time(1);
  flux180_data = flux180.data;
  nodes = 1:2;
  units_scale = 1e-4; % m^-2 > cm^-2 
  plot(hca,x_vec,abs(FVdv_edi_180)*units_scale,flux180_time+0.5*diff(flux180_time(1:2)),mean(flux180_data(:,nodes),2));  
  hca.YLabel.String = 'flux (cm^{-2}s^{-1})';
  legend(hca,'model','EDI','location','eastoutside')
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

axes_width = h(2).Position(3);

for ipanel = 1:npanels
  hca = h(ipanel);
  hca.Position(3) = axes_width;
  hca.Position(3);
  hca.XLim = x_vec([1 end]);
end
linkaxes(h,'x')

%%

% as the particle move into the potential field, the velocity changes

% potential energy in wave frame
u_es = @(phi) units.e*phi/units.eV;
U_POT = u_es(phi(k,X));
%Etot = E + U_POT;

% velocity grid
vmax = 3*vt;
nv = 500;
v = linspace(-vmax,vmax,nv);

[X,V] = meshgrid(x,v);
E_J = units.me*V.^2/2; % J
E = E_J/units.e; % res=Units.me*data.^2*10^6/2/Units.eV;
F0 = f(V,n,vt,vd,0);
%Fdrift = f(V,0.1*n,0.2*vt,20000e3,0);
Fphi = f(V,-0.2*n,vt,vph,phi(k,X));
Ftot = F0 + Fphi; % + Fdrift;

% energy and velocity of particle
% no wave
V_part = V; 
E_part = E;

% these are changed due to the wave

% plot
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % F without wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,F0);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Fphi);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
end

colormap(cn.cmap('blue_red'))

% see also
% open lic.phasespace
% open lic.phasespace_1

% Map phase space density, start with maxwellian, shift it to vph, add
% potential energy from wave, shift back to spacecraft frame, see how the
% maxwellian is modulated
%

units = irf_units;

% Plasma properties
n = [0.01]*1e6;
T = [400]; T_K = T*units.eV/units.kB; % use parallel temperature
vd = [0]; % m/s

% Physical parameters
vt = sqrt(2*units.e*T./units.me); % m/s
wp = sqrt(n.*units.e^2./(units.me*units.eps0)); % Hz
Le = units.c./wp;
Ld = vt./wp/sqrt(2);

% Wave properties
phi0 = 400; % V
vph = 8000e3; % m/s
E_vph = units.me*vph^2/2/units.eV;
wavelength = 8e3; % m
k = 2*pi/wavelength;
phi = @(k,x) phi0*sin(k*x); % same units as phi0
% x can be vph * t if we want to comapre to time series
% spatial grid (can also be time series, i.e. use directly)
t = 0:(1/8000/10):0.001;
x = t*vph;
nx = numel(x);

% 1D Maxwellian
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2);
f = @(v,n,vt,vd,phi) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2 - units.e*phi/(units.me*vt^2/2));

% potential energy in wave frame
u_es = @(phi) units.e*phi/units.eV;
U_POT = u_es(phi(k,X));
%Etot = E + U_POT;

% velocity grid
vmax = 3*vt;
nv = 500;
v = linspace(-vmax,vmax,nv);

[X,V] = meshgrid(x,v);
E_J = units.me*V.^2/2; % J
E = E_J/units.e; % res=Units.me*data.^2*10^6/2/Units.eV;
F0 = f(V,n,vt,vd,0);
%Fdrift = f(V,0.1*n,0.2*vt,20000e3,0);
Fphi = f(V,-0.2*n,vt,vph,phi(k,X));
Ftot = F0 + Fphi; % + Fdrift;

% energy and velocity of particle
% no wave
V_part = V; 
E_part = E;

% these are changed due to the wave

% plot
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % F without wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,F0);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Fphi);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
end

colormap(cn.cmap('blue_red'))

%% see also
% open lic.phasespace
% open lic.phasespace_1

% Map phase space density, start with maxwellian, shift it to vph, add
% potential energy from wave, shift back to spacecraft frame, see how the
% maxwellian is modulated
%

units = irf_units;

% Plasma properties
n = [0.01]*1e6;
T = [400]; T_K = T*units.eV/units.kB; % use parallel temperature
vd = [0]; % m/s

% Physical parameters
vt = sqrt(2*units.e*T./units.me); % m/s
wp = sqrt(n.*units.e^2./(units.me*units.eps0)); % Hz
Le = units.c./wp;
Ld = vt./wp/sqrt(2);

% Wave properties
phi0 = 400; % V
vph = 8000e3; % m/s
E_vph = units.me*vph^2/2/units.eV;
wavelength = 8e3; % m
k = 2*pi/wavelength;
phi = @(k,x) phi0*cos(k*x); % same units as phi0
% x can be vph * t if we want to comapre to time series
% spatial grid (can also be time series, i.e. use directly)
t = 0:(1/8000/10):0.001;
x = t*vph;
nx = numel(x);

% 1D Maxwellian
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2);
f = @(v,n,vt,vd,phi) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2 - units.e*phi/(units.me*vt^2/2));

% potential energy in wave frame
u_es = @(phi) units.e*phi/units.eV;
U_POT = u_es(phi(k,X));
%Etot = E + U_POT;

% velocity grid
vmax = 3*vt;
nv = 500;
v = linspace(-vmax,vmax,nv);

[X,V] = meshgrid(x,v);
E_J = units.me*V.^2/2; % J
E = E_J/units.e; % res=Units.me*data.^2*10^6/2/Units.eV;
F0 = f(V,n,vt,vd,0);
%Fdrift = f(V,0.1*n,0.2*vt,20000e3,0);
Fphi = f(V,-0.2*n,vt,vph,phi(k,X));
Ftot = F0 + Fphi; % + Fdrift;

% energy and velocity of particle
% no wave
V_part = V; 
E_part = E;

% these are changed due to the wave

% plot
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % F without wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,F0);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Fphi);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
end

colormap(cn.cmap('blue_red'))

% see also
% open lic.phasespace
% open lic.phasespace_1

% Map phase space density, start with maxwellian, shift it to vph, add
% potential energy from wave, shift back to spacecraft frame, see how the
% maxwellian is modulated
%

units = irf_units;

% Plasma properties
n = [0.01]*1e6;
T = [400]; T_K = T*units.eV/units.kB; % use parallel temperature
vd = [0]; % m/s

% Physical parameters
vt = sqrt(2*units.e*T./units.me); % m/s
wp = sqrt(n.*units.e^2./(units.me*units.eps0)); % Hz
Le = units.c./wp;
Ld = vt./wp/sqrt(2);

% Wave properties
phi0 = 400; % V
vph = 8000e3; % m/s
E_vph = units.me*vph^2/2/units.eV;
wavelength = 8e3; % m
k = 2*pi/wavelength;
phi = @(k,x) phi0*cos(k*x); % same units as phi0
% x can be vph * t if we want to comapre to time series
% spatial grid (can also be time series, i.e. use directly)
t = 0:(1/8000/10):0.001;
x = t*vph;
nx = numel(x);

% 1D Maxwellian
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2);
f = @(v,n,vt,vd,phi) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt^2 - units.e*phi/(units.me*vt^2/2));

% potential energy in wave frame
u_es = @(phi) units.e*phi/units.eV;
U_POT = u_es(phi(k,X));
%Etot = E + U_POT;

% velocity grid
vmax = 3*vt;
nv = 500;
v = linspace(-vmax,vmax,nv);

[X,V] = meshgrid(x,v);
E_J = units.me*V.^2/2; % J
E = E_J/units.e; % res=Units.me*data.^2*10^6/2/Units.eV;
F0 = f(V,n,vt,vd,0);
%Fdrift = f(V,0.1*n,0.2*vt,20000e3,0);
Fphi = f(V,-0.2*n,vt,vph,phi(k,X));
Ftot = F0 + Fphi; % + Fdrift;

% energy and velocity of particle
% no wave
V_part = V; 
E_part = E;

% these are changed due to the wave

% plot
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % F without wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,F0);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Fphi);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
end
if 1 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
end

colormap(cn.cmap('blue_red'))

%%
% plot
nrows = 5;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % E without wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,E);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E';
end
if 1 % E with wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,(V+vph)*1e-6,E);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E';
end
if 1 % U_POT
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,(V+vph)*1e-6,U_POT);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E';
end
if 1 % E-U_POT with wave
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,(V+vph)*1e-6,E-U_POT);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E';
end

for ipanel = 1:npanels
  hca = h(ipanel);
  hca.YLim = 70*[-1 1];
  %hca.CLim = 30000*[0 1]; 
  hca.CLim = (5*T+E_vph)*[0 1]; 
end
%%


u_ek_in_vph_frame = @(m,v,vp)(0.5*units.me*(v-vp).^2);

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
phi = @(k,x,o,t) 0.1*cos(k*x)+0*1;%.*(1+0.1*x);
e = -1;
m = 1;
t = 0;

V=V;

nPanels = 1;
for kk=1:nPanels; h(kk)=subplot(nPanels,1,kk); end
    
U_es = u_es(phi(k,X,o,t),e) + 0*u_es(phi(-k,X,o,t),e);
U_ek = u_ek_in_vph_frame(m,V,vp) + 0*u_ek_in_vph_frame(m,V,vp); 
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


