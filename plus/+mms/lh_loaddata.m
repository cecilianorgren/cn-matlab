% LH_LOADDATA Loads data.
%   LH_LOADDATA is a script that loads mms data needed to study lower 
%   hybrid waves with the matching method. 
%
%   The basic data needed is:
%       B - magnetic field
%       E - electric field
%       n - density
%
%   To compare with physical variables (length scales):
%       Te, Ti -> gyroradii, plasma beta
%       

tint = EpochTT(irf_time(toepoch([2015 08 15 12 30 00;2015 08 15 12 33 00])','epoch>utc'));
%tint = EpochTT(irf_time(toepoch([2015 09 15 12 30 00;2015 09 19 15 33 00])','epoch>utc'));
%tint = EpochTT(irf_time(toepoch([2015 08 15 13 00 01;2015 08 15 13 03 00])','epoch>utc'));
tint = EpochTT(irf_time(toepoch([2015 08 15 12 59 01;2015 08 15 13 03 00])','epoch>utc'));
tint = EpochTT(irf_time(toepoch([2015 08 15 12 59 01;2015 08 15 13 03 00])','epoch>utc'));

sc = 3;
% Magnetic field survey data
c_eval('[B?fg,dobj] = cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint(1));',sc)

% Magnetic field search coil data
% '20150828_1'
event = '20150815_1';
load_event_data;

% Electric field burst data
c_eval('[dslE?,dobj] = cn_get_ts(''mms?_edp_brst_ql_dce2d'',''mms3_edp_dce_xyz_dsl'',tint(1));',sc)

% Spacecraft potential
c_eval('[P?,dobjP?] = cn_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint(1));',sc);

% Ion skymap
c_eval('[disSkyMap?,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint(1));',sc)
c_eval('disPSD?.data = nanmean(nanmean(disSkyMap?.data,4),3);',sc)
c_eval('disPSD?.time = disSkyMap?.time;',sc);

% Electron skymap
c_eval('[desSkyMap?,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint(1));',sc)
c_eval('desPSD?.data = nanmean(nanmean(desSkyMap?.data,4),3);',sc)
c_eval('desPSD?.time = desSkyMap?.time;',sc);


%% Electron moments

c_eval('[~,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_des-moms'',[],tint(1));',sc)

c_eval('ne? = mms.variable2ts(get_variable(dobj,''mms?_des_numberDensity''));',sc)

c_eval('peXX = mms.variable2ts(get_variable(dobj,''mms?_des_PresXX'')); ',sc)
c_eval('peYY = mms.variable2ts(get_variable(dobj,''mms?_des_PresYY'')); ',sc)
c_eval('peZZ = mms.variable2ts(get_variable(dobj,''mms?_des_PresZZ'')); ',sc)
c_eval('pe? = irf.ts_scalar(peXX.time,(peXX.data + peYY.data + peZZ.data)/3);',sc)
c_eval('pe?.units = peXX.units;',sc)
c_eval('pe?.userData = peXX.userData;',sc)

c_eval('TeXX = mms.variable2ts(get_variable(dobj,''mms?_des_TempXX''));',sc)
c_eval('TeYY = mms.variable2ts(get_variable(dobj,''mms?_des_TempYY''));',sc)
c_eval('TeZZ = mms.variable2ts(get_variable(dobj,''mms?_des_TempZZ''));',sc)
c_eval('Te? = irf.ts_scalar(TeXX.time,(TeXX.data + TeYY.data + TeZZ.data)/3);',sc)
c_eval('Te?.units = TeXX.units;',sc)
c_eval('Te?.userData = TeXX.userData;',sc)

c_eval('veX = mms.variable2ts(get_variable(dobj,''mms?_des_bulkX''));',sc)
c_eval('veY = mms.variable2ts(get_variable(dobj,''mms?_des_bulkY''));',sc)
c_eval('veZ = mms.variable2ts(get_variable(dobj,''mms?_des_bulkZ''));',sc)
c_eval('ve? = irf.ts_vec_xyz(veX.time,[veX.data veY.data veZ.data]);',sc)
c_eval('ve?.units = veX.units;',sc)
c_eval('ve?.userData = veX.userData;',sc)

% Downsample electron moments
c_eval('fs = 1/(ne?.time(2)-ne?.time(1));',sc)
fny = fs/2;
c_eval('pe?_lowres = irf_filt(pe?,0,fny/2,fs,5);',sc)
c_eval('Te?_lowres = irf_filt(Te?,0,fny/2,fs,5);',sc)
c_eval('ne?_lowres = irf_filt(ne?,0,fny/2,fs,5);',sc)

%% Ion moments
c_eval('[~,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_dis-moms'',[],tint(1));',sc)

c_eval('ni? = mms.variable2ts(get_variable(dobj,''mms?_dis_numberDensity''));',sc)

c_eval('piXX = mms.variable2ts(get_variable(dobj,''mms?_dis_PresXX'')); ',sc)
c_eval('piYY = mms.variable2ts(get_variable(dobj,''mms?_dis_PresYY'')); ',sc)
c_eval('piZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_PresZZ'')); ',sc)
c_eval('pi? = irf.ts_scalar(piXX.time,(piXX.data + piYY.data + piZZ.data)/3);',sc)
c_eval('pi?.units = piXX.units;',sc)
c_eval('pi?.userData = piXX.userData;',sc)

c_eval('TiXX = mms.variable2ts(get_variable(dobj,''mms?_dis_TempXX''));',sc)
c_eval('TiYY = mms.variable2ts(get_variable(dobj,''mms?_dis_TempYY''));',sc)
c_eval('TiZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_TempZZ''));',sc)
c_eval('Ti? = irf.ts_scalar(TiXX.time,(TiXX.data + TiYY.data + TiZZ.data)/3);',sc)
c_eval('Ti?.units = TiXX.units;',sc)
c_eval('Ti?.userData = TiXX.userData;',sc)

c_eval('viX = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkX''));',sc)
c_eval('viY = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkY''));',sc)
c_eval('viZ = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkZ''));',sc)
c_eval('vi? = irf.ts_vec_xyz(viX.time,[viX.data viY.data viZ.data]);',sc)
c_eval('vi?.units = viX.units;',sc)
c_eval('vi?.userData = viX.userData;',sc)

% Downsample electron moments
c_eval('fs = 1/(ni?.time(2)-ni?.time(1));',sc)
fny = fs/2;
c_eval('pi?_lowres = irf_filt(pi?,0,fny/2,fs,5);',sc)
c_eval('Ti?_lowres = irf_filt(Ti?,0,fny/2,fs,5);',sc)
c_eval('ni?_lowres = irf_filt(ni?,0,fny/2,fs,5);',sc)
plotData = 0;
if plotData
    %%
    h = irf_plot(6);
    hca = irf_panel('Bfg');
    irf_plot(hca,B3fg)
    hca = irf_panel('Bsc');
    irf_plot(hca,B3sc)
    hca = irf_panel('E');
    irf_plot(hca,dslE3)
    hca = irf_panel('n');
    irf_plot(hca,ne3_lowres)
    hca = irf_panel('Electron PSD');
    specrec.p = desPSD3.data; specrec.t = desPSD3.time.epochUnix; specrec.f = 1:32; specrec.f_label = 'Energy level'; specrec.p_label = {'Electron phase','space density'};   
    irf_spectrogram(hca,specrec);
    hca = irf_panel('Ion PSD');
    specrec.p = disPSD3.data; specrec.t = disPSD3.time.epochUnix; specrec.f = 1:32; specrec.f_label = 'Energy level'; specrec.p_label = {'Ion phase','space density'};
    irf_spectrogram(hca,specrec);
    
    irf_zoom(h,'x',tint.epochUnix'+[-200 10])
    irf_zoom(h,'y')
    hca=gca;
    %irf_pl_mark(h,[tint(1).epochUnix tint(2).epochUnix])
    %irf_pl_mark(hca,[tint(1).epochUnix tint(2).epochUnix])
end