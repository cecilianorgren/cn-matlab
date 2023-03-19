%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 1;

localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
db_info = datastore('mms_db');   

tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% Time from time interval
tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');
tint_action = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');

% Time from file name
fileId = '20170725220853';

iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));

if 0
tint_all = irf.tint('2017-07-09T17:30:00.00Z/2017-07-09T17:35:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);
iFile = 1;
fileId = '20170709173053';
end

tint = [files(iFile-1).start files(iFile).stop] + [1 -1];

%tint = irf.tint('2017-07-09T17:30:00.00Z/2017-07-09T17:35:00.00Z');


% Event path
eventPath = ['/Users/' localuser '/Research/Events/mms_' fileId '/']; % for saving figures
eventPath = ['/Users/' localuser '/GoogleDrive/Research/Events/mms_' fileId '/']; % for saving figures
%matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms_' fileID '/'];
mkdir(eventPath)

%% Load data
% OMNI data, %% get omni pressure data
if 0
tint_omni = irf.tint('2017-07-25T20:00:00.00Z/2017-07-26T00:00:00.00Z');
tint_omni_long = irf.tint('2017-07-25T00:00:00.00Z/2017-07-26T00:00:00.00Z');
omni_solarwind = irf_get_data_omni(tint_omni.epochUnix,'n,v,P,Bz,Ms','omni_min'); % not working, dont know why
omni_index_min = irf_get_data_omni(tint_omni_long.epochUnix,'ae,al,au,omni_min'); % not working, dont know why
omni_index_hour = irf_get_data_omni(tint_omni_long.epochUnix,'ae,al,au,f10.7,dst,kp','omni_hour'); % not working, dont know why

ae = irf.ts_scalar(irf_time(omni_index_min(:,1),'epoch>EpochTT'),omni_index_min(:,2)); ae.name = 'AE';
al = irf.ts_scalar(irf_time(omni_index_min(:,1),'epoch>EpochTT'),omni_index_min(:,3)); al.name = 'AL';
au = irf.ts_scalar(irf_time(omni_index_min(:,1),'epoch>EpochTT'),omni_index_min(:,4)); au.name = 'AU';
end
%% Magnetic field
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',1:4);
% Electric field
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',1:4);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',1:4);
c_eval('gsmE? = c_coord_trans(''GSE'',''GSM'',gseE?);',1:4)
% Spacecraft potential
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',1:4);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',1:4);

%% HPCA
ic = 1:4;
c_eval('nHp?_brst = mms.get_data(''Nhplus_hpca_brst_l2'',tint,?);',ic);
c_eval('gsmVHp?_brst = mms.get_data(''Vhplus_gsm_hpca_brst_l2'',tint,?);',ic);
c_eval('gseVHp?_brst = c_coord_trans(''GSM'',''GSE'',gsmVHp?_brst);',ic)

%% Velocity
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)

% Pressure
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
% Pressure, FAC
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)
c_eval('Ti?perp = 0.5*(facTi?.yy+facTi?.zz);',ic)
c_eval('Te?perp = 0.5*(facTe?.yy+facTe?.zz);',ic)
c_eval('Ti?par = facTi?.xx;',ic)
c_eval('Te?par = facTe?.xx;',ic)
c_eval('Pi?perp = 0.5*(facPi?.yy+facPi?.zz);',ic)
c_eval('Pi?par = facPi?.xx;',ic)

% Spacecraft position
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)
%% Distributions
ic = 1;
% Distributions, FPI
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
% Remove all one-count "noise"
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_onecount = iPDist?; iPDist?_onecount.data = (iPDist?_onecount.data./iPDistErr?.data).^2;',ic)
%c_eval('iPDist?.data(iPDist?.data < iPDistErr?.data*1.01) = 0;',ic)

%c_eval('iPDist?_nobg = iPDist?;',ic)
%c_eval('iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)

c_eval('iPitch? = iPDist1_nobg.pitchangles(dmpaB?,12);',ic)
c_eval('iPitch?_nobg = iPDist1_nobg.pitchangles(dmpaB?,12);',ic)



tint_epd = [EpochTT('2017-07-25T22:04:44.822634526Z') EpochTT('2017-07-25T22:11:31.820198151Z')];
c_eval('eis_omni?_oplus = mms.get_data(''Omnifluxoxygen_epd_eis_brst_l2'',tint_epd,?);',1:4)
c_eval('eis_omni? = mms.get_data(''Omnifluxproton_epd_eis_brst_l2'',tint_epd,?);',1:4)
c_eval('eis_pa? = mms.get_data(''Pitchanglefluxproton_epd_eis_brst_l2'',tint_epd,?);',1:4)
c_eval('feeps_ion_omni? = mms.get_data(''Omnifluxion_epd_feeps_brst_l2'',tint_epd,?);',1:4)
c_eval('feeps_ion_pa? = mms.get_data(''Pitchanglefluxion_epd_feeps_brst_l2'',tint_epd,?);',1:4)

c_eval('feeps_pa? = mms.get_data(''Pitchanglefluxion_epd_feeps_brst_l2'',tint_epd,?);',1:4)
c_eval('feeps_ele_omni? = mms.get_data(''Omnifluxelectron_epd_feeps_brst_l2'',tint_epd,?);',1:4)


% Feeps ion omni, coverage not taken into account, 
omnidata1234 = (feeps_ion_omni1.data + ...
                feeps_ion_omni2.data(1:end-1,:) + ...
                feeps_ion_omni3.data + ...
                feeps_ion_omni4.data)/4;
feeps_ion_omni1234 =  PDist(feeps_ion_omni1.time,omnidata1234,'omni',feeps_ion_omni1.depend{1}); % energies keV -> eV            
feeps_ion_omni1234.siConversion = feeps_ion_omni1.siConversion;
feeps_ion_omni1234.units = feeps_ion_omni1.units;
feeps_ion_omni1234.species = feeps_ion_omni1.species;

% Feeps ion pitchangle
padata1234 = cat(4,feeps_ion_pa1.data,feeps_ion_pa2.data(1:end-1,:,:),feeps_ion_pa3.data,feeps_ion_pa4.data);  
padata1234 = nanmean(padata1234,4);
feeps_ion_pa1234 =  PDist(feeps_ion_pa1.time,padata1234,'pitchangle',feeps_ion_pa1.depend{1},feeps_ion_pa1.depend{2}); % energies keV -> eV            
feeps_ion_pa1234.siConversion = feeps_ion_pa1.siConversion;
feeps_ion_pa1234.units = feeps_ion_pa1.units;
feeps_ion_pa1234.species = feeps_ion_pa1.species;

% EIS, proton, mms1 has no data
omnidata1234 = (...
                eis_omni2.data(1:end-1,:) + ...
                eis_omni3.data + ...
                eis_omni4.data)/4;
eis_omni1234 =  PDist(eis_omni3.time,omnidata1234,'omni',eis_omni3.depend{1}); % energies keV -> eV            
eis_omni1234.siConversion = eis_omni3.siConversion;
eis_omni1234.units = eis_omni3.units;
eis_omni1234.species = eis_omni3.species;


padata1234 = cat(4,eis_pa2.data(1:end-1,:,:),eis_pa3.data,eis_pa4.data);
padata1234 = nanmean(padata1234,4);
eis_pa1234 =  PDist(eis_pa3.time,padata1234,'pitchangle',eis_pa3.depend{1},eis_pa3.depend{2}); % energies keV -> eV            
eis_pa1234.siConversion = eis_pa3.siConversion;
eis_pa1234.units = eis_pa3.units;
eis_pa1234.species = eis_pa3.species;

%% Derived quantities, par/perp, mag mom, beta, J, 
ic = 1:4;
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
% vExB
%c_eval('gsmVExB? = cross(gsmE?.resample(gsmB?.time),gsmB?)/gsmB?.abs/gsmB?.abs*1e3; gsmVExB?.units = '''';',ic) % km/s
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
% Magnetic field pressure
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
% Dynamic pressure
c_eval('PDi? = ne?.resample(gsmVi?)*1e6*0.5*units.mp*gsmVi?.x.^2*1e6*1e9; PDi?.name = ''Dynamic pressure''; PDi?.units =''nPa'';',ic)
% Plasma beta
c_eval('betae? = (gsePe?.trace.resample(gsePi?))/3/PB?.resample(gsePi?);',ic)
c_eval('betai? = (gsePi?.trace)/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = (gsePi?.trace+gsePe?.trace.resample(gsePi?))/3/PB?.resample(gsePi?); beta?.name = ''beta ie''; beta?.data(beta?.data<0) = NaN;',ic)
% Thermal speeds
c_eval('vti?perp = units.s*(1-1/(Ti?perp*units.e/(units.mp*units.c^2)+1).^2).^0.5*1e-3; ; vti?perp.units = ''km/s'';',ic)
c_eval('vte?perp = units.s*(1-1/(Te?perp*units.e/(units.me*units.c^2)+1).^2).^0.5*1e-3; ; vte?perp.units = ''km/s'';',ic)
c_eval('vti?par = units.s*(1-1/(Ti?par*units.e/(units.mp*units.c^2)+1).^2).^0.5*1e-3; ; vti?par.units = ''km/s'';',ic)
c_eval('vte?par = units.s*(1-1/(Te?par*units.e/(units.me*units.c^2)+1).^2).^0.5*1e-3; ; vte?par.units = ''km/s'';',ic)
% Magnetic moment based on thermal speeds
c_eval('mag_mome? = 0.5*units.me*vte?perp.^2*10^6/(gseB?.abs*1e-9)*1e9;  mag_mome?.units = ''nAm^2''; mag_mome?.name = ''magnetic moment'';',ic)
c_eval('mag_momi? = 0.5*units.me*vti?perp.^2*10^6/(gseB?.abs*1e-9)*1e9;  mag_momi?.units = ''nAm^2''; mag_momi?.name = ''magnetic moment'';',ic)
% Velocity, FAC
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)


c_eval('[gseVHp?par,gseVHp?perp] = irf_dec_parperp(gseB?,gseVHp?_brst); gseVHp?par.name = ''Vp par''; gseVHp?perp.name = ''Vp perp'';',ic)

% Currents from moments, use ne also for Ji 
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

%% Things requiring 4 sc.
c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?)); gseJxB?.name = ''JxB''; gseJxB?.units = ''nA/m^2 nT'';',1:4 )

c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.name = ''v_ixB'';',1:4)
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.name = ''v_exB'';',ic)

% Current from magnetic field
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');
% [JxB] = T A/m2


gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';

gseB.name = 'B'; gseDivB.units = 'nT';
gseDivPb.name = 'div P_B'; gseDivB.units = '...';
gseCurvB.name = 'curv B'; gseCurvB.units = '...';
gseJxB = JxB; gseJxB.name = 'JxB'; gseJxB.data = gseJxB.data; gseJxB.units = 'Am^-2 T';

c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)

gseBav = (gseB1.resample(gseB2) + gseB2.resample(gseB2) + gseB3.resample(gseB2) + gseB4.resample(gseB2))/4; gseBav.name = 'B1234'; 
gseEav = (gseE1.resample(gseE2) + gseE2.resample(gseE2) + gseE3.resample(gseE2) + gseE4.resample(gseE2))/4; gseEav.name = 'E1234'; % gseE1 not there?
[gseJcurlpar,gseJcurlperp] = irf_dec_parperp(gseBav,gseJcurl); gseJcurlpar.name = 'J curl par'; gseJcurlperp.name = 'J curl perp';


c_eval('gseR?brsttime = gseR?.resample(gseE?);',1:4)
[Ecurl,~,~,~,~,~] = c_4_j('gseR?brsttime','gseE?');
Ecurl.data = Ecurl.data;%*(units.mu0);
dtB = diff(gseBav.time-gseBav.time(1));
dBdt = diff(gseBav.data,1)./dtB;
dBdt = irf.ts_vec_xyz(gseBav.time(1:end-1)+0.5*dtB,dBdt);

% 4sc averages
Peav = (gsePe1.trace.resample(gsePe1) + gsePe2.trace.resample(gseVe1) + gsePe3.trace.resample(gseVe1) + gsePe4.trace.resample(gseVe1))/4/3; gsePeav.name = 'Pe1234';
Piav = (gsePi1.trace.resample(gsePi1) + gsePi2.trace.resample(gseVi1) + gsePi3.trace.resample(gseVi1) + gsePi4.trace.resample(gseVi1))/4/3; gsePiav.name = 'Pi1234';
PDiav = (PDi1.resample(PDi1) + PDi2.resample(PDi1) + PDi3.resample(PDi1) + PDi4.resample(PDi1))/4; PDiav.name = 'Pdyn_x_i1234';
PBav = (PB1.resample(PB1) + PB2.resample(PB1) + PB3.resample(PB1) + PB4.resample(PB1))/4; gsePiav.name = 'Pi1234';

gseVixBav = (gseVixB1.resample(gseVixB2) + gseVixB2.resample(gseVixB2) + gseVixB3.resample(gseVixB2) + gseVixB4.resample(gseVixB2))/4; gseVixBav.name = 'gseVixB1234';

%gseJxB = gseJcurl.cross(Bbrst); gseJxB.name = 'JxB'; gseJxB.units = 'nAm^-2 nT';

gseJav = (gseJ1.resample(gseJ2) + gseJ2.resample(gseJ2) + gseJ3.resample(gseJ2) + gseJ4.resample(gseJ2))/4; gseJav.name = 'J1234'; 
gseJiav = (gseJi1 + gseJi2.resample(gseJi1) + gseJi3.resample(gseJi1) + gseJi4.resample(gseJi1))/4;
gseJeav = (gseJe1 + gseJe2.resample(gseJe1) + gseJe3.resample(gseJe1) + gseJe4.resample(gseJe1))/4;

c_eval('JdotE? = gseJ?.dot(gseE?.resample(gseJi?))*1e3; JdotE?.name = ''JdotE''; JdotE?.units = ''nW/m^3'';',ic) % nA/m^2*mV/m = piko W

neav = (ne1.resample(ne1) + ne2.resample(ne1) + ne3.resample(ne1) + ne4.resample(ne1))/4; % gseE1 not there?
gseJxBne_mVm = (gseJxB)/(neav.resample(gseJxB)*1e6)/units.e*1e3; gseJxBne_mVm.name = 'JxB/ne';
gseJxBne_mVm.data(abs(gseJxBne_mVm.data)>100) = NaN;
%gseJxBne_mVm.data(abs(neav.data)<0.02,:) = NaN;

c_eval('gseJxBne?_mVm = (gseJxB?*1e-18)/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJxBne?_mVm.name = ''JxB/ne'';',ic)
%gseJxBne_mVm.data(abs(gseJxBne_mVm.data)>100) = NaN;


gseVexBav = (gseVexB1.resample(gseVexB2) + gseVexB2.resample(gseVexB2) + gseVexB3.resample(gseVexB2) + gseVexB4.resample(gseVexB2))/4; gseVexBav.name = 'gseVexB1234';

gseVeav = (gseVe1.resample(gseVe1) + gseVe2.resample(gseVe1) + gseVe3.resample(gseVe1) + gseVe4.resample(gseVe1))/4; gseVeav.name = 'gseVeB1234';
gseVeperpav = (gseVe1perp.resample(gseVe1) + gseVe2perp.resample(gseVe1) + gseVe3perp.resample(gseVe1) + gseVe4perp.resample(gseVe1))/4; gseVeav.name = 'gseVeperpB1234';
gseVeparav = (gseVe1par.resample(gseVe1par) + gseVe2par.resample(gseVe1) + gseVe3par.resample(gseVe1) + gseVe4par.resample(gseVe1))/4; gseVeav.name = 'gseVeB1234';

gseGradPe = divP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km'; gseGradPe.name = 'div Pe';
gseGradPi = divP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';
gseGradTe = divP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';
gseGradTi = divP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';
gseGradNe = c_4_grad('gseR?','ne?','grad');
gseGradNi = c_4_grad('gseR?','ni?','grad');

gseGradPene = gseGradPe/neav.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseGradPene.units = 'mV/m';
gseGradPene.data(abs(gseGradPene.data)>100) = NaN;

% Too noisy, I just keep it here
%gseGradPine = gseGradPi/neav.resample(gseGradPi.time)/units.e*1e-9*1e-6; gseGradPine.units = 'mV/m';
%gseGradPine.data(abs(gseGradPine.data)>100) = NaN;
%%

c_eval('iPitch? = iPDist1.pitchangles(dmpaB?,12);',ic)

disp('Done loading data.')

%% Reduced and pitchangle distributions
if 0
  %%
%c_eval('iPitch? = iPDist?.pitchangles(dmpaB?.resample(iPDist?),12);',ic)
ic = 1;
disp('Preparing reduced distributions.')
vint = [-Inf Inf];
elim = [200 40000];



c_eval('if1Dx?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[1 0 0],''vint'',vint);',ic)
c_eval('if1Dy?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dz?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[0 0 1],''vint'',vint);',ic)

%%
c_eval('if1Dx? = iPDist?.elim(elim).reduce(''1D'',[1 0 0],''vint'',vint);',ic)
c_eval('if1Dy? = iPDist?.elim(elim).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dz? = iPDist?.elim(elim).reduce(''1D'',[0 0 1],''vint'',vint);',ic)

%%
c_eval('if1Dx?_high = iPDist?.elim([10000 Inf]).reduce(''1D'',[1 0 0],''vint'',vint);',ic)
c_eval('if1Dy?_high = iPDist?.elim([10000 Inf]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dz?_high = iPDist?.elim([10000 Inf]).reduce(''1D'',[0 0 1],''vint'',vint);',ic)
c_eval('if1Dx?_low = iPDist?.elim([0 10000]).reduce(''1D'',[1 0 0],''vint'',vint);',ic)
c_eval('if1Dy?_low = iPDist?.elim([0 10000]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dz?_low = iPDist?.elim([0 10000]).reduce(''1D'',[0 0 1],''vint'',vint);',ic)


%%
vint = [-Inf Inf];
elim = [200 40000];
c_eval('if1DL? = iPDist?.elim(elim).reduce(''1D'',L,''vint'',vint);',ic)
c_eval('if1DM? = iPDist?.elim(elim).reduce(''1D'',M,''vint'',vint);',ic)
c_eval('if1DN? = iPDist?.elim(elim).reduce(''1D'',N,''vint'',vint);',ic)


c_eval('if1DL?_nobg = iPDist?_nobg.elim(elim).reduce(''1D'',L,''vint'',vint);',ic)
c_eval('if1DM?_nobg = iPDist?_nobg.elim(elim).reduce(''1D'',M,''vint'',vint);',ic)
c_eval('if1DN?_nobg = iPDist?_nobg.elim(elim).reduce(''1D'',N,''vint'',vint);',ic)

%%
%%
c_eval('if1Dy?_050 = iPDist?_nobg.elim([050 1e6]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dy?_100 = iPDist?_nobg.elim([100 1e6]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dy?_200 = iPDist?_nobg.elim([200 1e6]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)

c_eval('if1Dx?_nobg = iPDist?_nobg.reduce(''1D'',[1 0 0],''vint'',vint);',ic)
c_eval('if1Dy?_nobg = iPDist?_nobg.reduce(''1D'',[0 1 0],''vint'',vint);',ic)
c_eval('if1Dz?_nobg = iPDist?_nobg.reduce(''1D'',[0 0 1],''vint'',vint);',ic)

if 0
  %%
lowerelim = 100; % eV
c_eval('vPara = dmpaB?.resample(ePDist?).norm;',ic)
c_eval('ef1D?_par = ePDist?.elim([100 Inf]).tlim(tint_fered).reduce(''1D'',vPara,''vint'',vint,''scpot'',scPot?.resample(ePDist?));',ic)% reduced distribution along B
end
disp('Done preparing reduced distributions.')

%c_eval('gsmVe?.data(abs(gsmVe?.data)>20*1e3) = NaN;',ic)

%c_eval('[fitdata?,ts?] = funFitVDF(if1Dx?,''nPop'',2,''plot'',0,''guessprevious'',0,''X0'',[0.5e6 0 100e3 0.5e6 0 100e3]);',ic)
end
