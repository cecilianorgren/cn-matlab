tint = irf.tint('2015-10-16T10:32:50.00Z/2015-10-16T10:34:30.00Z'); % magnetosphere-magnetosheath-magnetosphere
ic = 1:4;

% xy -> yx (interchangeable)
% xz -> zx (interchangeable)
% yz -> zy (interchangeable)

%% Get pressure from fpi's des-moms
c_eval('Pexx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXX'',tint);',ic);
c_eval('Pexy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXY'',tint);',ic);
c_eval('Pexz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXZ'',tint);',ic);
c_eval('Peyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYY'',tint);',ic);
c_eval('Peyz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYZ'',tint);',ic);
c_eval('Pezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresZZ'',tint);',ic);

c_eval('Pdata? = nan(Pexx?brst.length,3,3);',ic)
c_eval(['Pdata?(:,1,1) = Pexx?brst.data;',...
        'Pdata?(:,1,2) = Pexy?brst.data;',...
        'Pdata?(:,1,3) = Pexz?brst.data;',...
        'Pdata?(:,2,1) = Pexy?brst.data;',...
        'Pdata?(:,2,2) = Peyy?brst.data;',...
        'Pdata?(:,2,3) = Peyz?brst.data;',...
        'Pdata?(:,3,1) = Pexz?brst.data;',...
        'Pdata?(:,3,2) = Peyz?brst.data;',...
        'Pdata?(:,3,3) = Pezz?brst.data;',...
        ],ic);

c_eval('P?fpi = irf.ts_tensor_xyz(Pexx?brst.time,Pdata?);',ic)

%% Get pressure from mms.psd_moments
c_eval('tmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l1b/des-dist/2015/10/16/mms?_fpi_brst_l1b_des-dist_20151016103254_v1.1.0.cdf'']);',ic);

c_eval('pdist? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_dist''));',ic);
c_eval('energy0? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy0'');',ic);
c_eval('energy1? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy1'');',ic);
c_eval('phi? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_phi''));',ic);
c_eval('theta? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_theta'');',ic);
c_eval('stepTable? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_stepTable_parity''));',ic);

c_eval('emoments? = mms.psd_moments(pdist?,phi?,theta?,stepTable?,energy0?,energy1?,P?brst,''electron'');',ic);

c_eval('Pet? = emoments?.P_psd.tlim(tint);',ic);

c_eval('Pdata? = nan(Pet?.length,3,3);',ic)
c_eval(['Pdata?(:,1,1) = Pet?.data(:,1);',...
        'Pdata?(:,1,2) = Pet?.data(:,2);',...
        'Pdata?(:,1,3) = Pet?.data(:,3);',...
        'Pdata?(:,2,1) = Pet?.data(:,2);',...
        'Pdata?(:,2,2) = Pet?.data(:,4);',...
        'Pdata?(:,2,3) = Pet?.data(:,5);',...
        'Pdata?(:,3,1) = Pet?.data(:,3);',...
        'Pdata?(:,3,2) = Pet?.data(:,5);',...
        'Pdata?(:,3,3) = Pet?.data(:,6);',...
        ],ic);
c_eval('P?irf = irf.ts_tensor_xyz(Pet?.time,Pdata?);',ic)

%% Calculate pressure gradient
ic = 1:4;
calPf1 = 1; calPf2 = 1; calPf3 = 1.08; calPf4 = 1;
c_eval('P? = {Pexx?brst*calPf?,Pexy?brst*calPf?,Pexz?brst*calPf?;Pexy?brst*calPf?,Peyy?brst*calPf?,Peyz?brst;Pexz?brst*calPf?,Peyz?brst*calPf?,Pezz?brst*calPf?};',ic)

gradP = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,P1,P2,P3,P4);
gradPe = gradP;
c_eval('tracePe? = (Pexx?brst + Peyy?brst + Pezz?brst)/3;')
avne = (ne1brst*calPf1 + ne2brst.resample(ne1brst.time)*calPf2 + ne3brst.resample(ne1brst.time)*calPf3 + ne4brst.resample(ne1brst.time)*calPf4)/4;
%gradP3 = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,P1,P2,P3,P4);


%%

c_eval('Texx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXX'',tint);',ic);
c_eval('Texy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXY'',tint);',ic);
c_eval('Texz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXZ'',tint);',ic);

%c_eval('Teyx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYX'',tint);',ic);
c_eval('Teyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYY'',tint);',ic);
c_eval('Teyz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYZ'',tint);',ic);

%c_eval('Tezx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZX'',tint);',ic);
%c_eval('Tezy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZY'',tint);',ic);
c_eval('Tezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZZ'',tint);',ic);
%%
if 0
  c_eval('Texx?filt = irf_filt(Texx?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
  c_eval('Texy?filt = irf_filt(Texy?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
  c_eval('Texz?filt = irf_filt(Texz?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
  c_eval('Teyy?filt = irf_filt(Teyy?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
  c_eval('Teyz?filt = irf_filt(Teyz?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
  c_eval('Tezz?filt = irf_filt(Tezz?brst,0,fsmoments?/4,fsmoments?,3);',[1:4]);
end
c_eval('T? = {Texx?brst,Texy?brst,Texz?brst;Texy?brst,Teyy?brst,Teyz?brst;Texz?brst,Teyz?brst,Tezz?brst};',ic)
gradT = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,T1,T2,T3,T4);
%% Construct pressure tensor
c_eval('Pdata? =  zeros(Pexx?brst.length,3,3); Pdata?(:,1:3,1) = [Pexx?brst.data, Pexy?brst.data, Pexz?brst.data]; Pdata?(:,1:3,2) = [Pexy?brst.data, Peyy?brst.data, Peyz?brst.data]; Pdata?(:,1:3,3) = [Pexz?brst.data, Peyz?brst.data, Pezz?brst.data];',ic)
c_eval('Pe?Ten = TSeries(Pexx?brst.time,Pdata?,''TensorOrder'',2,''repres'',{''x'',''y'',''z''},''repres'',{''x'',''y'',''z''});',ic)

%% Get pressure gradient
%%
ic = 3;
h = irf_plot(6);

hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('brst E');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[0.95 0.95]);

hca = irf_panel('P tensor');
c_eval('irf_plot(hca,{Pexx?brst,Pexy?brst,Pexz?brst,Peyy?brst,Peyz?brst,Pezz?brst},''comp'')',ic)
irf_legend(hca,{'P_{xx}','P_{xy}','P_{xz}','P_{yy}','P_{yz}','P_{zz}'},[0.95 0.95]);
hca.YLabel.String = {'P','(nPa)'};

hca = irf_panel('T tensor');
c_eval('irf_plot(hca,{Texx?brst,Texy?brst,Texz?brst,Teyy?brst,Teyz?brst,Tezz?brst},''comp'')',ic)
irf_legend(hca,{'T_{xx}','T_{xy}','T_{xz}','T_{yy}','T_{yz}','T_{zz}'},[0.95 0.95]);
hca.YLabel.String = {'T_{xx}','(eV)'};

hca = irf_panel('Pyy tensor 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{Peyy1brst,Peyy2brst,Peyy3brst,Peyy4brst},'comp')
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);
hca.YLabel.String = {'P_{yy}','(nPa)'};

hca = irf_panel('Pzz tensor 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{Pezz1brst,Pezz2brst,Pezz3brst,Pezz4brst},'comp')
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);
hca.YLabel.String = {'P_{zz}','(nPa)'};


%%
ic = 3;
h = irf_plot(3);

hca = irf_panel('P(xyz)x');
c_eval('irf_plot(hca,{Pexx?brst,Pexy?brst,Pexz?brst},''comp'')',ic)

hca = irf_panel('P(xyz)y');
c_eval('irf_plot(hca,{Pexy?brst,Peyy?brst,Peyz?brst},''comp'')',ic)

hca = irf_panel('P(xyz)z');
c_eval('irf_plot(hca,{Pexz?brst,Peyz?brst,Pezz?brst},''comp'')',ic)


%%
%c_eval('Pe?brst=irf.ts_vec_xyz(Pexx?brst.time,[Pexx?brst.data Peyy?brst.data Pezz?brst.data]);',ic)
%c_eval('Pe?Ten = TSeries(epoch,data4x3x3,''TensorOrder'',2,''repres'',{''x'',''y'',''z''},''repres'',{''x'',''y'',''z''})',ic)
