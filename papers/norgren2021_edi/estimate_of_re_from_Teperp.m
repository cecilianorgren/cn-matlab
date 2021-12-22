%% Define spacecraft and time intervals to be used
% Spacecraft id
ic = 1;

% Time interval of event
tint_burst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_burst = tint_burst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file

% Time intervals for modelling the distribution
tint_model = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); 
tint_model = irf.tint('2017-07-06T13:54:05.51Z/2017-07-06T13:54:05.63Z'); 

% Time interval for figure
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.620Z');

tint_warm = irf.tint('2017-07-06T13:54:05.00Z/2017-07-06T13:54:06.00Z');  
  
tint_cold = irf.tint('2017-07-06T13:54:13.000Z/2017-07-06T13:54:13.040Z');

%% Set up database
localuser = datastore('local','user');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS'); 
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

% MMS id and time interval
%mms_id = 1;
%tint = irf.tint('2017-07-06T13:54:05.520Z/2017-07-06T13:54:05.640Z');
%c_eval('eDist = ePDist?;',mms_id)
%c_eval('eDist_bgremoved = eDist_nobg?.tlim(tint);',mms_id)

% Load EH properties
% observed/measured properties (konrad)
% data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
% obs_eh_properties = data_tmp.EH_properties;
% obs_lpp = obs_eh_properties.Lpp; % peak to peak length
% obs_potential = obs_eh_properties.potential;
% obs_vtrap_all = sqrt(2*units.e*obs_potential/units.me);
% obs_potential_max = obs_eh_properties.max_potential;
% obs_velocity = obs_eh_properties.vel;
% obs_neh = numel(obs_velocity);

% EDI parameters
% E_edi = 500; % eV
% v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
% dE_edi = 25; % eV
% 
% E_edi_plus = E_edi + dE_edi;
% E_edi_minus = E_edi - dE_edi;
% v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
% v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
% v_edi_plusminus = v_edi_plus-v_edi_minus;
% dv_edi_minus = v_edi_minus - v_edi;
% dv_edi_plus = v_edi_plus - v_edi;
% dv_edi = dv_edi_plus - dv_edi_minus; % m/s

% Remove background electrons
nPhoto = 00;
nSecond = 20;
if 0
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint, ...
             'Nphotoe_art', nPhoto, 'nSecondary', nSecond, 'ZeroNaN', 0);
elseif 0
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint,...
             'ZeroNaN', 0);  
end

%% Make reduced distribution
strTint = [irf_time(tint(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vint = [-Inf Inf];
vg = (-40:1:40)*1e3;
c_eval('eDist = ePDist?.tlim(tint_figure);',ic)
c_eval('eDist = eDist_nobg?.tlim(tint_figure);',ic)

           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 40;
nMC = 500;
tic; ef1D_par        = eDist.reduce('1D',ePara,        'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp1      = eDist.reduce('1D',ePerp1,       'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp2      = eDist.reduce('1D',ePerp2,       'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1   = eDist.reduce('2D',ePara,ePerp1, 'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2   = eDist.reduce('2D',ePara,ePerp2, 'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

% Cold populations
c_eval('eDist = ePDist?.tlim(tint_cold);',ic)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 40;
nMC = 500;
tic; ef1D_par_cold        = eDist.reduce('1D',ePara,        'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp1_cold      = eDist.reduce('1D',ePerp1,       'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp2_cold      = eDist.reduce('1D',ePerp2,       'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1_cold   = eDist.reduce('2D',ePara,ePerp1, 'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2_cold   = eDist.reduce('2D',ePara,ePerp2, 'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2_cold = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc


% tic; ef1D_bgremoved = eDist_bgremoved.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
% tic; ef2D_parperp1_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
% tic; ef2D_parperp2_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
% tic; ef2D_perp1perp2_bgremoved = eDist_bgremoved.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

% Make pitch angle spectrograms
% ePitch = eDist.pitchangles(dmpaB1.resample(eDist),12);
% ePitch_bgremoved = eDist_bgremoved.pitchangles(dmpaB1.resample(eDist_bgremoved),12); 

%% Plot distributions and fit for ESW time
units = irf_units;
f = @(v,n,vd,vt) n*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
fk = @(v,n,vd,vt,k) f(v,n,vd,vt).*(1+(v-vd).^2./k/vt.^2).^(-k-1);
v = ef1D_par.depend{1}(1,:)*1e3;

if 0 % not used?
% warm, parallel
ntot = 0.034*1e6;
n = ntot;
T = 200; % warm populations is 200 eV
vd = 0e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fwarm = f(v,n,vd,vt);

% hot, parallel
ntot = 0.006*1e6;
n = ntot;
T = 2000; % warm populations is 200 eV
vd = 2000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fhot = f(v,n,vd,vt);
end

% warm, parallel
ntot = 0.034*1e6;
n = ntot;
T = 400; % warm populations is 200 eV
vd = -5000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fwarm = f(v,n,vd,vt);

% hot, for parallel
ntot = 0.006*1e6;
n = ntot;
T = 1300; % warm populations is 200 eV
vd = 10000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fhot = f(v,n,vd,vt);

% warm, for perp1
ntot = 0.034*1e6;
n = ntot;
T = 300; % warm populations is 200 eV
vd = 000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fwarm_perp1 = f(v,n,vd,vt);

ntot = 0.034*1e6;
nk = ntot*1.9;
Tk = 1000; % warm populations is 200 eV
vdk = 000e3;
vtk = sqrt(2*units.e*Tk./units.me); % m/s
kappa = 0.3;
fwarm_perp1_k = fk(v,nk,vdk,vtk,kappa);

% warm, for perp2
ntot = 0.034*1e6;
n = ntot;
T = 250; % warm populations is 200 eV
vd = 000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fwarm_perp2 = f(v,n,vd,vt);

ntot = 0.034*1e6;
nk = ntot*1.9;
Tk = 800; % warm populations is 200 eV
vdk = 000e3;
vtk = sqrt(2*units.e*Tk./units.me); % m/s
kappa = 0.3;
fwarm_perp2_k = fk(v,nk,vdk,vtk,kappa);

% hot, perp
ntot = 0.006*1e6;
n = ntot;
T = 1300; % warm populations is 200 eV
vd = 0000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par.depend{1}(1,:)*1e3;
fhot_perp = f(v,n,vd,vt);

nrows = 2;
ncols = 3;
isub = 1;
hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_par.depend{1}(1,:)*1e-3,mean(ef1D_par.data,1),...
  ef1D_par.depend{1}(1,:)*1e-3,fwarm,...
  ef1D_par.depend{1}(1,:)*1e-3,fhot,...
  ef1D_par.depend{1}(1,:)*1e-3,fwarm+fhot)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_perp1.depend{1}(1,:)*1e-3,mean(ef1D_perp1.data,1),...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp1,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fhot_perp,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp1+fhot_perp,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp1_k,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp1_k+fhot_perp)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_perp2.depend{1}(1,:)*1e-3,mean(ef1D_perp2.data,1),...
  ef1D_perp2.depend{1}(1,:)*1e-3,fwarm_perp2,...
  ef1D_perp2.depend{1}(1,:)*1e-3,fhot_perp,...
  ef1D_perp2.depend{1}(1,:)*1e-3,fwarm_perp2+fhot_perp,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp2_k,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp2_k+fhot_perp)


if 1 % f*v^2
hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_par.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_par.data,1),...
  v,v.^2.*fwarm,...
  v,v.^2.*fwarm_k)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_perp1.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_perp1.data,1),...
  v,v.^2.*fwarm_perp1,...
  v,v.^2.*fwarm_perp1_k)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_perp2.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_perp2.data,1),...
  v,v.^2.*fwarm_perp2,...
  v,v.^2.*fwarm_perp2_k)
end

c_eval('h(?).YLim = [0 0.4e-2];',1:3)
c_eval('h(?).XLim = [-40 40];',1:3)
legend({'obs','fit warm','fit hot'})

c_eval('h(?).YScale = ''log'';',1:3)
c_eval('h(?).YLim = [1e-5 0.4e-2];',1:3)

%% Plot distributions and fit for lobe time
units = irf_units;
f = @(v,n,vd,vt) n*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
v = ef1D_par_cold.depend{1}(1,:)*1e3;

if 0 % not used?
% warm, parallel
ntot = 0.034*1e6;
n = ntot;
T = 200; % warm populations is 200 eV
vd = 0e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par_cold.depend{1}(1,:)*1e3;
fwarm = f(v,n,vd,vt);

% hot, parallel
ntot = 0.006*1e6;
n = ntot;
T = 2000; % warm populations is 200 eV
vd = 2000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par_cold.depend{1}(1,:)*1e3;
fhot = f(v,n,vd,vt);
end

% warm, parallel
ntot = 0.04*1e6;
n = ntot;
T = 70; % warm populations is 200 eV
vd = -1000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par_cold.depend{1}(1,:)*1e3;
fwarm = f(v,n,vd,vt);
% kappa
ntot = 0.04*1e6;
nk = ntot*2;
Tk = 200; % warm populations is 200 eV
vdk = -1000e3;
vtk = sqrt(2*units.e*Tk./units.me); % m/s
kappa = 0.2;
fwarm_k = fk(v,nk,vdk,vtk,kappa);


% warm, for perp1
ntot = 0.04*1e6;
n = ntot;
T = 50; % warm populations is 200 eV
vd = 000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par_cold.depend{1}(1,:)*1e3;
fwarm_perp1 = f(v,n,vd,vt);
% kappa
ntot = 0.04*1e6;
nk = ntot*2;
Tk = 200; % warm populations is 200 eV
vdk = 000e3;
vtk = sqrt(2*units.e*Tk./units.me); % m/s
kappa = 0.3;
fwarm_perp1_k = fk(v,nk,vdk,vtk,kappa);

% warm, for perp2
ntot = 0.04*1e6;
n = ntot;
T = 70; % warm populations is 200 eV
vd = 000e3;
vt = sqrt(2*units.e*T./units.me); % m/s
v = ef1D_par_cold.depend{1}(1,:)*1e3;
fwarm_perp2 = f(v,n,vd,vt);
% kappa
ntot = 0.04*1e6;
nk = ntot*2.2;
Tk = 200; % warm populations is 200 eV
vdk = 000e3;
vtk = sqrt(2*units.e*Tk./units.me); % m/s
kappa = 0.3;
fwarm_perp2_k = fk(v,nk,vdk,vtk,kappa);


nrows = 2;
ncols = 3;
isub = 1;

hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_par_cold.depend{1}(1,:)*1e-3,mean(ef1D_par_cold.data,1),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,fwarm,...
  ef1D_par_cold.depend{1}(1,:)*1e-3,fwarm_k)
hca.XLabel.String = 'v (10^3 km/s)';
hca.YLabel.String = 'f_e (s/m^4)';


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_perp1_cold.depend{1}(1,:)*1e-3,mean(ef1D_perp1_cold.data,1),...
  ef1D_perp1_cold.depend{1}(1,:)*1e-3,fwarm_perp1,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp1_k)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_perp2_cold.depend{1}(1,:)*1e-3,mean(ef1D_perp2_cold.data,1),...
  ef1D_perp2_cold.depend{1}(1,:)*1e-3,fwarm_perp2,...
  ef1D_perp1.depend{1}(1,:)*1e-3,fwarm_perp2_k)

if 1 % f*v^2
hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_par_cold.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_par_cold.data,1),...
  v,v.^2.*fwarm,...
  v,v.^2.*fwarm_k)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_perp1_cold.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_perp1_cold.data,1),...
  v,v.^2.*fwarm_perp1,...
  v,v.^2.*fwarm_perp1_k)


hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
v = ef1D_perp2_cold.depend{1}(1,:)*1e-3;
plot(hca,...
  v,v.^2.*mean(ef1D_perp2_cold.data,1),...
  v,v.^2.*fwarm_perp2,...
  v,v.^2.*fwarm_perp2_k)
end
if 0 % cumsum of f
hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(mean(ef1D_par_cold.data,1)),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm_k))

hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(mean(ef1D_perp1_cold.data,1)),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm_perp1),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm_perp1_k))

hca = subplot(nrows,ncols,isub); h(isub) = hca; isub = isub + 1;
plot(hca,...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(mean(ef1D_perp2_cold.data,1)),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm_perp2),...
  ef1D_par_cold.depend{1}(1,:)*1e-3,cumsum(fwarm_perp2_k))
end

c_eval('h(?).YLim = [0 0.6e-2];',1:3)
c_eval('h(?).XLim = [-30 30];',1:3)
legend({'obs','fit max','fit kappa'},'location','best','box','off')

c_eval('h(?).YScale = ''log'';',1:3)
c_eval('h(?).YLim = [1e-6 1e-2];',1:3)

%% Calculate Tperp in parallel speed range around ESWs
tint = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
strTint = [irf_time(tint(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vlimpar = -8500 + 5000*[-1 1]; % km/s
%vlimpar = 1000000*[-1 1]; % km/s
vint = [-Inf Inf];
%vg = (-100:2:100)*1e3;
vgpar = vlimpar; % km/s
vgperp = (-50:1:50)*1e3; % km/s
c_eval('eDist = ePDist?.tlim(tint);',ic)
eDist = eDist(1);
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 00;
nMC = 500;
orient = [ePara.data; ePerp1.data; ePerp2.data];
%ePDistCart = eDist.elim([lowerelim inf]).rebin('cart',{vgpar,vgperp,vgperp},orient,'scpot',scpot);
ePDistCart = eDist.elim([lowerelim inf]).rebin('cart',{vgpar,vgperp,vgperp},orient,'scpot',scpot*0);

% Calculate temperature within given velocity interval
units = irf_units;
m = units.me;
f = squeeze(ePDistCart.data); % s^3cm^-6
v_scale = 1e3/1e-2; % km/s - > cm/s
dVX = diff(ePDistCart.ancillary.vx_edges)*v_scale;
[VY,VZ] = ndgrid(ePDistCart.depend{2}*v_scale,ePDistCart.depend{3}*v_scale);
[dVY,dVZ] = ndgrid(diff(ePDistCart.ancillary.vy_edges)*v_scale,diff(ePDistCart.ancillary.vz_edges)*v_scale);
mom_n = sum(sum(f.*dVX.*dVY.*dVZ))
mom_vy = sum(sum(f.*VY.*dVX.*dVY.*dVZ))
mom_vz = sum(sum(f.*VZ.*dVX.*dVY.*dVZ))
v_perp = sqrt((VY-mom_vy).^2 + (VZ-mom_vz).^2)*1e-2; % cm/s -> m/s
T_perp = sum(sum(m*f.*v_perp.^2.*dVX.*dVY.*dVZ));
%Tmatdv = Tmat.*dVY.*dVZ;
%sumT = sum(Tmatdv(:));
%%
% Plot results
hca = subplot(1,1,1);
pcolor(hca,ePDistCart.depend{2}*1e-3,ePDistCart.depend{3}*1e-3,log10(squeeze(ePDistCart.data))')
colormap(hca,pic_colors('candy4'))
shading(hca,'flat')
hca.XLabel.String = 'v_{perp,1} (10^3 km/s)';
hca.YLabel.String = 'v_{perp,2} (10^3 km/s)';
hca.Title.String = sprintf('%.0f<v_{||}<%.0f km/s',vlimpar(1),vlimpar(2));
hcb = colorbar('peer',hca);
hcb.YLabel.String = sprintf('log_{10} f (%s)',ePDistCart.units);