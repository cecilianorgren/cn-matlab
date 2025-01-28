tint = irf.tint('2015-12-02T01:14:14/2015-12-02T01:15:14Z');
%tint = irf.tint('2015-10-16T10:33:10.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere +- 10 s
ic = 1:4;
%% Electric field
disp('Loading electric field...')
%Eoffs1 = 2.05; Eoffs2 = 2.58; Eoffs3 = 2.53; Eoffs4 = 1.40;
c_eval('tic; dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
%c_eval('dslE?brst.data(:,1)=dslE?brst.data(:,1)-Eoffs?;',ic)

%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l1b_hmfe'',''mms?_edp_hmfe_sensor'',tint); toc',ic);
%c_eval('dslE?3Dbrst=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

%% Magnetic field
disp('Loading magnetic field...')
tic
c_eval('tic; dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint); toc',1:4);

%%
tref = irf.tint('2015-10-16T10:33:18.00Z',0.1);
c_eval('Bref? = dmpaB?brst.tlim(tref).data(1,:);') 
c_eval('dmpaB?brstRemOff = dmpaB?brst+Bref1-Bref?;') % Remove offset
%c_eval('dmpaB?brst = dmpaB?brstRemoff;') % Remove offset
%% 
disp('Loading spacecraft potential...')
c_eval('P?brst=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint);',ic);

% Magnetic field (searchcoil)
%c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_064900.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);
%c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_065200.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);
%%
% Electron moments
disp('Loading electron moments...')
c_eval('ne?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_numberDensity'',tint);',ic);

c_eval('vex?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkX'',tint);',ic);
c_eval('vey?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkY'',tint);',ic);
c_eval('vez?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkZ'',tint);',ic);
c_eval('ve?brst=irf.ts_vec_xyz(vex?brst.time,[vex?brst.data vey?brst.data vez?brst.data]);',ic)
c_eval('ve?brst.name = ''mms? ve brst'';')

c_eval('Pexx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXX'',tint);',ic);
c_eval('Pexy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXY'',tint);',ic);
c_eval('Pexz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXZ'',tint);',ic);

c_eval('Peyx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYX'',tint);',ic);
c_eval('Peyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYY'',tint);',ic);
c_eval('Peyz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYZ'',tint);',ic);

c_eval('Pezx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresZX'',tint);',ic);
c_eval('Pezy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresZY'',tint);',ic);
c_eval('Pezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresZZ'',tint);',ic);

c_eval('Pe?brst=irf.ts_vec_xyz(Pexx?brst.time,[Pexx?brst.data Peyy?brst.data Pezz?brst.data]);',ic)
%c_eval('Pe?Ten = TSeries(epoch,data4x3x3,''TensorOrder'',2,''repres'',{''x'',''y'',''z''},''repres'',{''x'',''y'',''z''})',ic)
%c_eval('Pdata? =  zeros(Pexx?brst.length,3,3); Pdata?(:,1:3,1) = [Pexx?brst.data, Pexy?brst.data, Pexz?brst.data]; Pdata?(:,1:3,2) = [Peyx?brst.data, Peyy?brst.data, Peyz?brst.data]; Pdata?(:,1:3,3) = [Pezx?brst.data, Pezy?brst.data, Pezz?brst.data];',ic)
%c_eval('Pdata? =  []; Pdata?(:,1:3,1) = [Pexx?brst.data, Pexx?brst.data, Pexx?brst.data]; Pdata?(:,1:3,2) = [Peyy?brst.data, Peyy?brst.data, Peyy?brst.data]; Pdata?(:,1:3,3) = [Pezz?brst.data, Pezz?brst.data, Pezz?brst.data];',ic)
%c_eval('Pe?Ten = TSeries(Pexx?brst.time,Pdata?,''TensorOrder'',2,''repres'',{''x'',''y'',''z''},''repres'',{''x'',''y'',''z''})',ic)

%c_eval('Texx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXX'',tint);',ic);
%c_eval('Teyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYY'',tint);',ic);
%c_eval('Tezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZZ'',tint);',ic);

c_eval('Texx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXX'',tint);',ic);
c_eval('Texy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXY'',tint);',ic);
c_eval('Texz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXZ'',tint);',ic);

c_eval('Teyx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYX'',tint);',ic);
c_eval('Teyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYY'',tint);',ic);
c_eval('Teyz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYZ'',tint);',ic);

c_eval('Tezx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZX'',tint);',ic);
c_eval('Tezy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZY'',tint);',ic);
c_eval('Tezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZZ'',tint);',ic);

c_eval('Te?brst=irf.ts_vec_xyz(Texx?brst.time,[Texx?brst.data Teyy?brst.data Tezz?brst.data]);',ic)

c_eval('fse = 1/(ne?brst.time(2)-ne?brst.time(1));',ic); fny = fse/2;
c_eval('Pe?_lowres = irf_filt(Pe?brst,0,fny/2,fse,5);',ic)
c_eval('Te?_lowres = irf_filt(Te?brst,0,fny/2,fse,5);',ic)
c_eval('ne?_lowres = irf_filt(ne?brst,0,fny/2,fse,5);',ic)

c_eval('[PeXXp?,PeXYp?,PeXZp?,PeYYp?,PeYZp?,PeZZp?] = mms.rotate_tensor_fac(Pexx?brst,Pexy?brst,Pexz?brst,Peyy?brst,Peyz?brst,Pezz?brst,dmpaB?brstRemOff);',ic)
c_eval('[TeXXp?,TeXYp?,TeXZp?,TeYYp?,TeYZp?,TeZZp?] = mms.rotate_tensor_fac(Texx?brst,Texy?brst,Texz?brst,Teyy?brst,Teyz?brst,Tezz?brst,dmpaB?brstRemOff);',ic)

c_eval('Te?par = TeZZp?;',1:4); c_eval('Tepar?_lowres = irf_filt(Te?par,0,fny/2,fse,5);',ic) 
c_eval('Te?perp = (TeXXp?+TeYYp?)/2;',1:4); c_eval('Teperp?_lowres = irf_filt(Te?perp,0,fny/2,fse,5);',ic) 


% Ion moments
%%
disp('Loading ion moments...')
c_eval('ni?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',tint);',ic);

c_eval('vix?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkX'',tint);',ic);
c_eval('viy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkY'',tint);',ic);
c_eval('viz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkZ'',tint);',ic);
c_eval('vi?brst=irf.ts_vec_xyz(vix?brst.time,[vix?brst.data viy?brst.data viz?brst.data]);',ic)

c_eval('Pixx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresXX'',tint);',ic);
c_eval('Piyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresYY'',tint);',ic);
c_eval('Pizz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresZZ'',tint);',ic);
c_eval('Pi?brst=irf.ts_vec_xyz(Pixx?brst.time,[Pixx?brst.data Piyy?brst.data Pizz?brst.data]);',ic)

c_eval('Tixx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempXX'',tint);',ic);
c_eval('Tiyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempYY'',tint);',ic);
c_eval('Tizz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempZZ'',tint);',ic);
c_eval('Ti?brst=irf.ts_vec_xyz(Tixx?brst.time,[Tixx?brst.data Tiyy?brst.data Tizz?brst.data]);',ic)

c_eval('fsi = 1/(ni?brst.time(2)-ni?brst.time(1));',ic); fny = fsi/2;
c_eval('Pi?_lowres = irf_filt(Pi?brst,0,fny/2,fsi,5);',ic)
c_eval('Ti?_lowres = irf_filt(Ti?brst,0,fny/2,fsi,5);',ic)
c_eval('ni?_lowres = irf_filt(ni?brst,0,fny/2,fsi,5);',ic)
c_eval('vi?_lowres = irf_filt(vi?brst,0,fny/4,fsi,5);')
%% c_eval('ni?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',tint);',ic);
tinttmp = tint;
%tint = irf.tint('2015-10-16T10:28:00.00Z/2015-10-16T10:40:00.00Z');
disp('Loading electron distribution data...')
c_eval('tic; desDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint); toc',[1 2 3 4]);
disp('Loading ion distribution data...')
c_eval('tic; disDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint); toc',[1 2 3 4]);
%%
if 0
dbRootPath = datastore('mms_db','local_file_db_root');
c_eval('tmpDobj = dataobj([dbRootPath ''/mms?/fpi/brst/l1b/dis-dist/2015/10/16/mms?_fpi_brst_l1b_dis-dist_20151016103000_v1.0.0.cdf'']); disDist? = mms.variable2ts(get_variable(tmpDobj,''mms?_dis_brstSkyMap_dist'')); disDist? = disDist?.tlim(tint);')
c_eval('tmpDobj = dataobj([dbRootPath ''/mms?/fpi/brst/l1b/des-dist/2015/10/16/mms?_fpi_brst_l1b_des-dist_20151016103000_v1.0.0.cdf'']); desDist? = mms.variable2ts(get_variable(tmpDobj,''mms?_des_brstSkyMap_dist'')); desDist? = desDist?.tlim(tint);')
%c_eval('desDist? = mms.variable2ts(mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint));',ic);
end
%%
%tmpDataObj = dataobj('/Volumes/SAMSUNG/data/mms1/fpi/brst/l1b/des-dist/2015/10/16/mms1_fpi_brst_l1b_des-dist_20151016103000_v0.2.0.cdf');
%dist = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_brstSkyMap_dist'));
%disttemp = dist.tlim(tint);

%% Construct variables
disp('Calculating ExB velocity, speeds, currents, frequencies, length scales...')
% Current from 4sc magnetic field: assuming GSE and DMPA are the same coordinate system.
c_eval('gseR?brsttime = gseR?.resample(dmpaB?brst);',1:4)
[jbrst,divBbrst,Bbrst,jxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','dmpaB?brst');
jbrst = irf.ts_vec_xyz(jbrst.time,jbrst.data);
jbrst.data = jbrst.data*1e9; jbrst.units = 'nAm^{-2}';
jbrst.time = EpochTT(jbrst.time); jbrst.name = '4sc current density';
%%
% Currents from moments
e = 1.6022e-19;
c_eval('je? = -e*ne?brst*ve?brst*1e3*1e6*1e9; je?.units = ''nA/m^2'';',ic);
c_eval('ji? = e*ni?brst*vi?brst*1e3*1e6*1e9; ji?.units = ''nA/m^2'';',ic);
c_eval('ji?a = irf.ts_vec_xyz(ji?.time(1:2:end,:),ji?.data(1:2:end,:));',ic);
c_eval('ji?b = irf.ts_vec_xyz(ji?.time(2:2:end,:),ji?.data(2:2:end,:));',ic);
c_eval('ji? = (ji?a+ji?b.resample(ji?a.time))*0.5; ji?.units = ''nA/m^2'';',ic);
c_eval('jtot? = (je?+ji?.resample(je?.time));',ic);

%%
% Gradient of electron pressure
% start of this time interval is a relatively quiet period
tref = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z')+[1 0]; 
c_eval('Pref? = Pe?_lowres.tlim(tref).data(1,:);') % Pa
c_eval('Pe?rel1 = Pe?_lowres+Pref1-Pref?;') % Remove offset
c_eval('r? = [Pe1rel1.time.epochUnix gseR?.resample(Pe1rel1.time).data Pe?rel1.abs.resample(Pe1rel1.time).data/3];')
for ii = 1:size(r1,1)
  gradPhie(ii,1)=r1(ii,1);
  gradPhie(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
gradPe = irf.ts_vec_xyz(irf_time(gradPhie(:,1),'epoch>utc'),gradPhie(:,2:4));
gradPe.units = 'Pa/km';
gradPe.name = 'grad(Pe)';

%%
% Gradient of electron pressure
% start of this time interval is a relatively quiet period
tref = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z')+[1 0]; 
c_eval('Pe?_ = Pe?_lowres;',[1 2 4])
Pe3_ = Pe3_lowres*1.1;%*1.1; % calibrated to mms3 to time where they observe similar field
c_eval('r? = [Pe1_.time.epochUnix gseR?.resample(Pe1_.time).data Pe?_.x.resample(Pe1_.time).data Pe?_.y.resample(Pe1_.time).data Pe?_.z.resample(Pe1_.time).data];')
divPeTmp = [];
divPeTmpVal = [];
for ii = 1:size(r1,1)
  gradPeTmp(ii,1)=r1(ii,1);
  gradPeTmpVal = c_4_gradvec(r1,r2,r3,r4);
  gradPeTmp(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
%[divPe,curlPe] = c_4_gradvec(r1,r2,r3,r4);
gradPe = irf.ts_vec_xyz(irf_time(gradPeTmp(:,1),'epoch>utc'),gradPeTmp(:,2:4));
gradPe.units = 'Pa/km';
gradPe.name = 'grad(Pe)';
%% Ions
tref = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z')+[1 0]; 
c_eval('Pi?_ = Pi?_lowres;',[1 2 3 4])
%Pi3_ = Pi3_lowres*1.1;%*1.1; % calibrated to mms3 to time where theyobserve similar field
c_eval('r? = [Pi1_.time.epochUnix gseR?.resample(Pi1_.time).data Pi?_.x.resample(Pi1_.time).data Pi?_.y.resample(Pi1_.time).data Pi?_.z.resample(Pi1_.time).data];')
gradPiTmp = [];
gradPiTmpVal = [];
for ii = 1:size(r1,1)
  gradPiTmp(ii,1)=r1(ii,1);
  gradPiTmpVal = c_4_gradvec(r1,r2,r3,r4);
  gradPiTmp(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
%[divPi,curlPi] = c_4_gradvec(r1,r2,r3,r4);
gradPi = irf.ts_vec_xyz(irf_time(gradPiTmp(:,1),'epoch>utc'),gradPiTmp(:,2:4));
gradPi.units = 'Pa/km';
gradPi.name = 'grad(Pi)';
%%
% Gradient of ion pressure, NOT GOOD
% start of this time interval is a relatively quiet period
tref = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z')+[1 0]; 
c_eval('Piref? = Pi?_lowres.tlim(tref).data(1,:);') % Pa
c_eval('Pi?rel1 = Pi?_lowres+Piref1-Piref?;') % Remove offset
Pi3rel1 = Pi3rel1*0.82; % Correct form
Pi4rel1 = Pi4rel1*1.1; % Correct form
c_eval('r? = [Pi1rel1.time.epochUnix gseR?.resample(Pi1rel1.time).data Pi?rel1.abs.resample(Pi1rel1.time).data/3];')
for ii = 1:size(r1,1)
  gradPhii(ii,1)=r1(ii,1);
  gradPhii(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
gradPi = irf.ts_vec_xyz(irf_time(gradPhii(:,1),'epoch>utc'),gradPhii(:,2:4));
gradPi.units = 'Pa/km';
gradPi.name = 'grad(Pi)';

%%
% Different representations of the electric field, for momentum balance
% equation and Ohm's law
e = 1.6022e-19;
avNe = (ne1brst + ne2brst.resample(ne1brst.time) + ne3brst.resample(ne1brst.time) + ne4brst.resample(ne1brst.time))/4;
%%
avB = (dmpaB1brst.resample(dmpaB1brst.time) + dmpaB2brst.resample(dmpaB1brst.time) + dmpaB3brst.resample(dmpaB1brst) + dmpaB4brst.resample(dmpaB1brst))/4;
c_eval('dmpaB?brst.units = ''nT'';',1:4);
avBbrst = (dmpaB1brst + dmpaB2brst.resample(dmpaB1brst.time) + dmpaB3brst.resample(dmpaB1brst) + dmpaB4brst.resample(dmpaB1brst))/4;
javxBav = j.cross(avB.resample(j.time));
javxBavbrst = j.cross(avBbrst.resample(j.time));
c_eval('jxB? = jtot?.cross(dmpaB?brst.resample(jtot?.time));')
jxBavmom = (jxB1 + jxB2.resample(jxB1) + jxB3.resample(jxB1) + jxB4.resample(jxB1))/4;
%%
avNe = (ne1_lowres + ne2_lowres.resample(ne1_lowres.time) + ne3_lowres.resample(ne1_lowres.time)*1.1 + ne4_lowres.resample(ne1_lowres.time))/4;
gradPene = gradPe*1e-8/avNe/e*1e-6; gradPene.units = 'mV/m';
c_eval('gradPene? = gradPe*1e-8/ne?_lowres/e*1e-6; gradPene?.units = ''mV/m'';')
avNi = (ni1_lowres + ni2_lowres.resample(ni1_lowres.time) + ni3_lowres.resample(ni1_lowres.time)*1.1 + ni4_lowres.resample(ni1_lowres.time))/4;
gradPini = gradPi*1e-8/avNi/e*1e-6; gradPini.units = 'mV/m';
c_eval('vexB? = ve?brst.resample(dmpaB?brst).cross(dmpaB?brst);')
c_eval('vexB?mVm = vexB?*1e-3; vexB?.units = ''mV/m'';')
c_eval('vixB?mVm = vixB?*1e-3; vixB?.units = ''mV/m'';')

c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')

%%

ic = 1:4;
% ExB velocity
c_eval('vExB?brst = cross(dslE?brst,dmpaB?brst.resample(dslE?brst.time))/dmpaB?brst.abs.resample(dslE?brst.time)/dmpaB?brst.abs.resample(dslE?brst.time)*1e3;',ic) % km/s

% Speeds
%c_eval('Bmat = double(dmpaB?brst.abs.data); nemat = ne?_lowres.resample(dmpaB?brst.time).data; Temat = Te?_lowres.abs.resample(dmpaB?brst.time).data; Timat = Ti?_lowres.abs.resample(dmpaB?brst.time).data;',ic)
c_eval('vte?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Vte''); vte?brst = irf.ts_scalar(dmpaB?brst.time,vte?brst);',ic)
c_eval('vtp?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Vtp''); vtp?brst = irf.ts_scalar(dmpaB?brst.time,vtp?brst);',ic)
c_eval('vA?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Va''); vA?brst = irf.ts_scalar(dmpaB?brst.time,vA?brst);',ic)

% Frequencies
c_eval('flh?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Flh''); flh?brst = irf.ts_scalar(dmpaB?brst.time,flh?brst);',ic)
c_eval('fce?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Fce''); fce?brst = irf.ts_scalar(dmpaB?brst.time,fce?brst);',ic)
c_eval('fcp?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Fcp''); fcp?brst = irf.ts_scalar(dmpaB?brst.time,fcp?brst);',ic)
c_eval('fpe?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Fpe''); fpe?brst = irf.ts_scalar(dmpaB?brst.time,fpe?brst);',ic)
c_eval('fpp?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Fpp''); fpp?brst = irf.ts_scalar(dmpaB?brst.time,fpp?brst);',ic)

% Length scales
c_eval('Lp?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Li''); Lp?brst = irf.ts_scalar(dmpaB?brst.time,Lp?brst);',ic)
c_eval('Le?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Le''); Le?brst = irf.ts_scalar(dmpaB?brst.time,Le?brst);',ic)
c_eval('Ld?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Ld''); Ld?brst = irf.ts_scalar(dmpaB?brst.time,Ld?brst);',ic)
c_eval('re?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Roe''); re?brst = irf.ts_scalar(dmpaB?brst.time,re?brst);',ic)
c_eval('rp?brst = irf_plasma_calc(dmpaB?brst.abs.data,ne?_lowres.resample(dmpaB?brst.time).data,0,Te?_lowres.abs.resample(dmpaB?brst.time).data,Ti?_lowres.abs.resample(dmpaB?brst.time).data,''Rop''); rp?brst = irf.ts_scalar(dmpaB?brst.time,rp?brst);',ic)

%% Wave spectrogram
