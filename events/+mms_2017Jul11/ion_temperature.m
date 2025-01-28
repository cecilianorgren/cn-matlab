%% Load data
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');
ic = 1;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); % torbert
%tint = irf.tint('2017-08-04T10:00:03.00Z/2017-08-04T10:04:03.00Z'); % hakon, shorter
%tint = irf.tint('2017-08-04T10:00:03.00Z/2017-08-04T10:01:03.00Z'); % hakon
%tint = irf.tint('2018-07-05T20:21:03.000Z/2018-07-05T20:25:52.000Z'); % hakon 2
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
% c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('hplusOmni? = mms.get_data(''Omnifluxhplus_hpca_brst_l2'',tint,?);',ic)
c_eval('oplusOmni? = mms.get_data(''Omnifluxoplus_hpca_brst_l2'',tint,?);',ic)
c_eval('Thplus? = mms.get_data(''Thplus_dbcs_hpca_brst_l2'',tint,?);',ic)
c_eval('nhplus? = mms.get_data(''Nhplus_hpca_brst_l2'',tint,?);',ic)
c_eval('noplus? = mms.get_data(''Noplus_hpca_brst_l2'',tint,?);',ic)
c_eval('gseTi?part = mms.get_data(''partTi_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('ni?part = mms.get_data(''partNi_fpi_brst_l2'',tint,?);',ic)
c_eval('gseTis?part = irf.ts_scalar(gseTi?part.time,squeeze(gseTi?part.data(:,1,:,1)+gseTi?part.data(:,2,:,2)+gseTi?part.data(:,3,:,3))/3);',ic)
c_eval('Ti?part = PDist(gseTi?part.time,squeeze(gseTi?part.data(:,1,:,1)+gseTi?part.data(:,2,:,2)+gseTi?part.data(:,3,:,3))/3,''moms-tens0'',gseTi?part.userData.energy.data);',ic)
c_eval('Ti?partdiff = PDist(gseTi?part.time,diff([zeros(gseTi?part.length,1) squeeze(gseTi?part.data(:,1,:,1)+gseTi?part.data(:,2,:,2)+gseTi?part.data(:,3,:,3))/3],1,2),''moms-tens0'',gseTi?part.userData.energy.data(:,2:end));',ic)
c_eval('eisi? = mms.get_data(''Omnifluxion_epd_eis_brst_l2'',tint,?);',ic)
c_eval('eisi?fast = mms.get_data(''Omnifluxion_epd_eis_fast_l2'',tint,?);',ic)
c_eval('feepsi? = mms.get_data(''Omnifluxion_epd_feeps_brst_l2'',tint,?);',ic)

%% Loading feeps data
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)

%%   Example to plot FEEPS/e data received from FEEPS team (Drew Turner, drew.lawson.turner@gmail.com). 
%   - Data:
%       1. OMNI particle flux of FEEPS-electron: "*_flux_omni_avg.txt";
%       2. energy table [15; keV]: "*_flux_omni_avg_v.txt";
%       3. pitch angle distribution of one energy level (e1 - e15): "*_flux_pad_e?_smth.txt";
%       4. pitch angle [19; -5 - 185 deg]: "*_flux_pad_e?_smth_v.txt";
%   - History:
%       1. wyli, 2016-03-11, irfu.
%       2. wyli, 2016-10-27, irfu; update code, IDL/tplot2txt.pro 

%% FEEPS
%filepath = '/Volumes/Fountain/Data/MMS/mms1/feeps/brst/l2/electron/2018/07/05/mms4_feeps_brst_l2_electron_20180705202323_v6.1.3.cdf';
filepath = '/Volumes/Fountain/Data/MMS/mms1/feeps/brst/l2/electron/2018/07/05/mms1_feeps_brst_l2_electron_20180705202323_v6.1.3.cdf';
dobj = dataobj(filepath);
time = get_variable(dobj,'epoch'); time = EpochTT(time.data);

pa = get_variable(dobj,'mms1_epd_feeps_brst_l2_electron_pitch_angle'); pa = pa.data; % [-180,180]
spinsectnum = get_variable(dobj,'mms1_epd_feeps_brst_l2_electron_spinsectnum'); % [0 63]

electron_energy = get_variable(dobj,'electron_energy');
electron_energy_lower_bound = get_variable(dobj,'electron_energy_lower_bound');
electron_energy_upper_bound = get_variable(dobj,'electron_energy_upper_bound');
%% EDP-EIS
dobj = dataobj('/Volumes/Fountain/Data/MMS/mms1/epd-eis/brst/l2/electronenergy/2017/07/11/mms1_epd-eis_brst_l2_electronenergy_20170711222923_v4.0.3.cdf');
dobj = dataobj('/Volumes/Fountain/Data/MMS/mms1/epd-eis/brst/l2/extof/2018/07/03/mms1_epd-eis_brst_l2_extof_20180703151343_v4.0.4.cdf');
dobj = dataobj('/Volumes/Fountain/Data/MMS/mms1/epd-eis/brst/l2/phxtof/2018/07/03/mms1_epd-eis_brst_l2_phxtof_20180703151343_v4.0.4.cdf');
% phxtof - pulse-hight vs (x) time-of-flight
% higher pulse hight suggests heavier ions
% extof - energy vs (x) time-of-flight
%% HPCA
ic = 1;
c_eval('hpcai_obj? = dataobj(''/Volumes/Fountain/Data/MMS/mms?/hpca/brst/l2/ion/2018/07/05/mms?_hpca_brst_l2_ion_20180705195423_v4.1.1.cdf'');', ic);
c_eval('HpfluxV? = get_variable(hpcai_obj?,''mms?_hpca_hplus_flux'');', ic);
c_eval('HpfluxT? = get_ts(hpcai_obj?,''mms?_hpca_hplus_flux'');', ic);
c_eval('HepfluxV? = get_variable(hpcai_obj?,''mms?_hpca_heplus_flux'');', ic);
c_eval('HepfluxT? = get_ts(hpcai_obj?,''mms?_hpca_heplus_flux'');', ic);
c_eval('HeppfluxV? = get_variable(hpcai_obj?,''mms?_hpca_heplusplus_flux'');', ic);
c_eval('HeppfluxT? = get_ts(hpcai_obj?,''mms?_hpca_heplusplus_flux'');', ic);    
c_eval('OpfluxV? = get_variable(hpcai_obj?,''mms?_hpca_oplus_flux'');', ic);
c_eval('OpfluxT? = get_ts(hpcai_obj?,''mms?_hpca_oplus_flux'');', ic);
c_eval('energyHp = HpfluxV?.DEPEND_2.data;', ic);
c_eval('energyHep = HepfluxV?.DEPEND_2.data;', ic);
c_eval('energyHepp = HeppfluxV?.DEPEND_2.data;', ic);
c_eval('energyOp = OpfluxV?.DEPEND_2.data;', ic);

% 3.0. fix mean0 error.
c_eval('HpfluxT?.data(HpfluxT?.data<=0) = NaN;', ic);    
c_eval('HepfluxT?.data(HepfluxT?.data<=0) = NaN;', ic);
c_eval('HeppfluxT?.data(HeppfluxT?.data<=0) = NaN;', ic);
c_eval('OpfluxT?.data(OpfluxT?.data<=0) = NaN;', ic);
%       3.1. HpfluxT1
c_eval('specHp=struct(''t'',HpfluxT?.time.epochUnix);', ic);
c_eval('specHp.p = squeeze(irf.nanmean(HpfluxT?.data, 2));', ic);
specHp.p_label={'flux','[1/cc s sr eV]'};
specHp.f_label={'Energy', '[eV]'};
specHp.f = single(energyHp);
%       3.2. HepfluxT1
c_eval('specHep=struct(''t'',HepfluxT?.time.epochUnix);', ic);
c_eval('specHep.p = squeeze(irf.nanmean(HepfluxT?.data, 2));', ic);
specHep.p_label={'flux','[1/cc s sr eV]'};
specHep.f_label={'Energy', '[eV]'};
specHep.f = single(energyHep);
%       3.3. HeppfluxT1
c_eval('specHepp=struct(''t'',HeppfluxT?.time.epochUnix);', ic);
c_eval('specHepp.p = squeeze(irf.nanmean(HeppfluxT?.data, 2));', ic);
specHepp.p_label={'flux','[1/cc s sr eV]'};
specHepp.f_label={'Energy', '[eV]'};
specHepp.f = single(energyHepp);
%       3.4. OpfluxT1
c_eval('specOp=struct(''t'',OpfluxT?.time.epochUnix);', ic);
c_eval('specOp.p = squeeze(irf.nanmean(OpfluxT?.data, 2));', ic);
specOp.p_label={'flux','[1/cc s sr eV]'};
specOp.f_label={'Energy', '[eV]'};
specOp.f = single(energyOp);    
specHepp.f = single(energyHepp);

%% HPCA

HpfluxV = get_variable(dobj,'mms1_hpca_hplus_flux');
HpfluxT = get_ts(dobj,'mms1_hpca_heplus_flux');
specrec = 1;

vars = dobj.vars;
for ivar = 1:numel(vars)
  hpca.(vars{ivar}) = get_variable(dobj,vars{ivar});
  
end

time = EpochTT(hpca.Epoch.data) + 1e-3*double(0.5*(hpca.Epoch_PLUS.data + hpca.Epoch_MINUS.data));
energy = hpca.mms1_hpca_ion_energy.data;
psd = hpca.mms1_hpca_hplus_phase_space_density.data;
azimuth = hpca.mms1_hpca_azimuth_angles_degrees.data;
polar = hpca.mms1_hpca_polar_anode_number.data;

%%
    ic = 1;
    Tint = irf.tint('2015-11-01T03:37:00.00Z/2015-11-01T03:41:10.00Z'); 
    feeps_dr = '/Users/wyli/wyli/work/data/feeps/20151101/feeps_ascii/';
    mmsx_omni_avg_fn = 'mmsx_feeps_brst_flux_omni_avg.txt';
    mmsx_omni_energy_fn = 'mmsx_feeps_brst_flux_omni_avg_v.txt';
    mmsx_omni_fn = 'mmsx_feeps_brst_flux_omni.txt';    
    c_eval('mmsx_pad_e?_fn = ''mmsx_feeps_brst_flux_pad_e?_smth.txt'';', 1: 15);
    c_eval('mmsx_pad_e?_PA_fn = ''mmsx_feeps_brst_flux_pad_e?_smth_v.txt'';', 1: 15);
    
%%  2. load MMSX-FEEPS data
    % 2.1. omni        
    formatomni=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 15])];
    IDomni = fopen([feeps_dr mmsx_omni_fn], 'r');
    IDomni_avg = fopen([feeps_dr mmsx_omni_avg_fn], 'r');    
    tmpData_omni = textscan(IDomni, formatomni);
    tmpData_omni_avg = textscan(IDomni_avg, formatomni);
    fclose(IDomni);         fclose(IDomni_avg);
    % 2.2. omni energy
    energy = load([feeps_dr mmsx_omni_energy_fn]);
    % 2.3. PAD E1-15
    formatpad=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 19])];
    c_eval('IDpad_e? = fopen([feeps_dr mmsx_pad_e?_fn], ''r'');', 1: 15);
    c_eval('tmpData_pad_e? = textscan(IDpad_e?, formatpad);', 1: 15);
    c_eval('fclose(IDpad_e?);', 1: 15);
    % 2.4. PAD E1-15 pitch angle
    c_eval('pa_e? = load([feeps_dr mmsx_pad_e?_PA_fn]);', 1: 15);

%%  3. load Bxyz
    c_eval('gseB=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);
    
%%  4. make
    % 4.1. omni & omni_avg data
    tomni = irf_time([tmpData_omni{1, 1}, tmpData_omni{1, 2}, tmpData_omni{1, 3}, ...
        tmpData_omni{1, 4}, tmpData_omni{1, 5}, tmpData_omni{1, 6}], 'vector6>epochtt');
    Data_omni = [tmpData_omni{1, 7: 21}];
    tomni_avg = irf_time([tmpData_omni_avg{1, 1}, tmpData_omni_avg{1, 2}, tmpData_omni_avg{1, 3}, ...
        tmpData_omni_avg{1, 4}, tmpData_omni_avg{1, 5}, tmpData_omni_avg{1, 6}], 'vector6>epochtt');
    Data_omni_avg = [tmpData_omni_avg{1, 7: 21}];    
    % 4.2. pad data
    c_eval('tpad_e? = irf_time([tmpData_pad_e?{1, 1}, tmpData_pad_e?{1, 2}, tmpData_pad_e?{1, 3}, tmpData_pad_e?{1, 4}, tmpData_pad_e?{1, 5}, tmpData_pad_e?{1, 6}], ''vector6>epochtt'');', 1: 15);
    c_eval('Data_pad_e? = [tmpData_pad_e?{1, 7: 25}];', 1: 15);
    % 4.3. make omni spectrum structure;
    Data_omni((Data_omni < 0.)) = NaN;
    specomni =struct('t', tomni.epochUnix);
    specomni.p = double(Data_omni); 
    specomni.p_label={'dF','#/(cm^2 s sr keV)'};
    specomni.f_label={''};
    specomni.f = single(energy);
    % 4.4. make omni_avg spectrum structure;  
    Data_omni_avg((Data_omni_avg < 0.)) = NaN;    
    specomni_avg =struct('t', tomni_avg.epochUnix);
    specomni_avg.p = double(Data_omni_avg); 
    specomni_avg.p_label={'dF','#/(cm^2 s sr keV)'};
    specomni_avg.f_label={''};
    specomni_avg.f = single(energy);
    % 4.5. make pitch angle spectrum structure;  
    c_eval('Data_pad_e?(Data_pad_e?<0) = NaN;', 1: 15);
    c_eval('specPAD_e? = struct(''t'', tpad_e?.epochUnix);', 1: 15);
    c_eval('specPAD_e?.p = double(Data_pad_e?);', 1: 15);
    c_eval('specPAD_e?.p_label={''dF'',''#/(cm^2 s sr keV)''};', 1: 15);
    c_eval('specPAD_e?.f_label={''''};', 1: 15);
    c_eval('specPAD_e?.f = single(pa_e?);', 1: 15);
    % 4.6. 52 keV flux
    dflux52 = irf.ts_scalar(tomni, Data_omni(:, 2));
   
%%  5. plot
    h = irf_plot(6,'newfigure');
    xSize=900; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.15;
    c_eval('set(h(?), ''position'', [0.12 0.965-? * ywidth xwidth ywidth]);', 1: 6);
    
    % 5.1. Bxyz
    h(1) = irf_panel('Bgse');
    irf_plot(h(1), gseB, 'LineWidth', 1.5);
    ylabel(h(1),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{x} ',' B_{y} ',' B_{z}'},[0.88 0.15])
    irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    
    % 5.2. omni;
    h(2) = irf_panel('omni');
    irf_spectrogram(h(2), specomni, 'log');
    caxis(h(2),[1 5]);
    ylabel(h(2), {'W_{e}', '[keV]'},'Interpreter','tex');
    irf_legend(h(2),'(b)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(2).YScale = 'log';
    irf_zoom(h(2), 'y', [45 400]);    
    
    % 5.3. omni;
    h(3) = irf_panel('omni_avg');
    irf_spectrogram(h(3), specomni_avg, 'log');
    caxis(h(3), [0 1.8]);
    ylabel(h(3), {'W_{e}^{avg}', '[keV]'},'Interpreter','tex');
    irf_legend(h(3),'(c)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    h(3).YScale = 'log';
    irf_zoom(h(3), 'y', [45 400]);
    
    % 5.4. PAD - e2;
    h(4) = irf_panel('PAD_e2');
    irf_spectrogram(h(4), specPAD_e2, 'log');
    caxis(h(4),[0 2]);
    ylabel(h(4), {'52 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(4),'(d)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    h(4).YTick = [0 45 90 135 180];
    
    % 5.5. PAD - e3;
    h(5) = irf_panel('PAD_e3');    
    irf_spectrogram(h(5), specPAD_e3, 'log');
    caxis(h(5),[0 1.5]);
    ylabel(h(5), {'70 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(5),'(e)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(5).YTick = [0 45 90 135 180];

    % 5.6. PAD - e5;
    h(6) = irf_panel('PAD_e5');    
    irf_spectrogram(h(6), specPAD_e5, 'log');
    caxis(h(6),[0 1.2]);
    ylabel(h(6), {'107 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(6),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(6).YTick = [0 45 90 135 180];    
    
    % 5.X. global set
    load('caa/cmap.mat');
    colormap(h(2),cmap);
    colormap(h(3),cmap); 
    colormap(h(4),cmap); 
    colormap(h(5),cmap); 
    colormap(h(6),cmap);     
    irf_plot_axis_align(h(1: 6));         % align the width of all panels
    irf_zoom(h(1: 6),'x',Tint);       % zoom all panels to the same time interval
    tinthl1 = irf_time('2015-11-01T03:37:30.0000Z', 'utc>epoch');    
    tinthl2 = irf_time('2015-11-01T03:40:50.8000Z', 'utc>epoch');
    irf_pl_mark(h(1:6), tinthl1, 'red', 'LineStyle', '--', 'LineWidth', 2.5);    
    irf_pl_mark(h(1:6), tinthl2, 'blue', 'LineStyle', '--', 'LineWidth', 2.5);
    title(h(1), 'MMSX FEEPS/e');
%%

% IDL/tplot2txt.pro: convert FEEPS-tplot data to *.txt data; need SPEDAS toolbox
%   PRO FEEPS_20151101
%   filenames = 'C:\Users\wyli\Desktop\20151101\mmsx_data_4Wenya_01Nov2015.tplot'
%   cd, 'C:\Users\wyli\Desktop\20151101\feeps_ascii'        ; windows system
%   tplot_restore ,filenames=filenames, all=all, sort=sort, get_tvars=get_tvars
%   tplot_ascii, 'mmsx_feeps_brst_flux_omni_avg'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_omni'   
%   tplot_ascii, 'mmsx_fgm_Bcomps_gsm'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e1'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e2'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e3'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e4'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e5'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e6'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e7'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e8'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e9'              ; ... 10 - 15  
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e1_smth'       ; 1-15
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e2_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e3_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e4_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e5_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e6_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e7_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e8_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e9_smth'        ; ... 10 - 15
%   tplot_ascii, 'mmsx_feeps_brst_elec_pas'
%   end

%% Figure
npanels = 6;
%h = irf_plot(npanels);
[h,h2] = initialize_combined_plot(npanels,2,1,0.6,'vertical');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z},'comp');
  hca.YLabel.String = {'B (GSE)','(nT)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % iPDist deflux omni  
  hca = irf_panel('i DEF omni');  
  irf_spectrogram(hca,iPDist1.deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 1 % Ti part
  %%
  figure;
  
  hca = irf_panel('i Ti part');  
  irf_spectrogram(hca,Ti1partdiff.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  Tsum = irf.ts_scalar(Ti1part.time,cumsum(Ti1part.data,2));
  irf_plot(hca,{gseTi1.trace/3,Tsum},'comp');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 1 % iPDist dpflux omni feeps 
  hca = irf_panel('i PEF omni feeps');  
  irf_spectrogram(hca,feepsi1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 1 % iPDist dpflux omni EIS
  hca = irf_panel('i PEF EIS omni');  
  irf_spectrogram(hca,Omnifluxion_epd_eis_brst_l2.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};  
end
if 1 % iPDist dpflux omni  FPI
  hca = irf_panel('i PEF omni');  
  irf_spectrogram(hca,iPDist1.dpflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i^{FPI}','(eV)'};  
end
if 1 % iPDist dpflux omni, HPCA  
  hca = irf_panel('i DPF omni HPCA');
  irf_spectrogram(hca,hplusOmni1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,Thplus1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_{H+}^{HPCA}','(eV)'}; 
  irf_legend(hca,{'HCPA'},[0.02,0.02],'color',[0 0 0])
end
if 0 % oPDist deflux omni, HPCA  
  hca = irf_panel('o DPF omni HPCA');  
  irf_spectrogram(hca,oplusOmni1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,Toplus1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_{O+}','(eV)'}; 
  irf_legend(hca,{'HCPA'},[0.02,0.02],'color',[0 0 0])
end
if 0 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  irf_spectrogram(hca,ePDist1.deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTe1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 0 % Ti, Te
  hca = irf_panel('T');
  set(hca,'ColorOrder',mms_colors('2134'))
  irf_plot(hca,{gseTi1.trace/3,Thplus1.trace/3,gseTe1.trace/3},'comp');
  hca.YLabel.String = {'T','(eV)'};
  irf_legend(hca,{'fpi i','hpca h+','e'},[0.98 0.98])
end
if 0 % ni,ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('2134'))
  irf_plot(hca,{ni1,nhplus1,noplus1,ne1},'comp');
  hca.YLabel.String = {'n','(cm^{-3})'};
  irf_legend(hca,{'fpi i','hpca h+','hpca o+','fpi e'},[0.98 0.98])
end
if 0 % gseE
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.x,gseE1.y,gseE1.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 0 % gseE, lowpass
  hca = irf_panel('E lowpass');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{plotE.x,plotE.y,plotE.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 0 % gseVixB
  hca = irf_panel('vixB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVixB1.x,-gseVixB1.y,-gseVixB1.z},'comp');
  hca.YLabel.String = {'-v_i\times B (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 0 % gseE, lowpass, vvixB x
  hca = irf_panel('E vixb x');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.x,plotE.x,-gseVixB1.x},'comp');
  hca.YLabel.String = {'E_x (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_x','E_x','-v_ixB_x'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 0 % gseE, lowpass, vvixB y
  hca = irf_panel('E vixb y');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.y,plotE.y,-gseVixB1.y},'comp');
  hca.YLabel.String = {'E_y (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_y','E_y','-v_ixB_y'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 0 % gseE, lowpass, vvixB z
  hca = irf_panel('E vixb z');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.z,plotE.z,-gseVixB1.z},'comp');
  hca.YLabel.String = {'E_z (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_z','E_z','-v_ixB_z'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end

irf_plot_axis_align
irf_zoom(h,'x',tint)

hca = h2(1);
omni_fpi.e = iPDist1.depend{1}(1,:);
omni_fpi.d = nanmean(iPDist1.dpflux.omni.data,1);
omni_hpca.e = hplusOmni1.depend{1}(:);
omni_hpca.d = nanmean(hplusOmni1.data,1);
omni_eis.e = Omnifluxion_epd_eis_brst_l2.depend{1}(1,:);
omni_eis.d = nanmean(Omnifluxion_epd_eis_brst_l2.data,1)*1e-3; % 1/keV -> /eV
omni_feeps.e = feepsi1.depend{1}(:);
omni_feeps.d = nanmean(feepsi1.data,1)*1e-3; % 1/keV -> /eV
loglog(hca,omni_fpi.e,omni_fpi.d,omni_hpca.e,omni_hpca.d,omni_eis.e,omni_eis.d,omni_feeps.e,omni_feeps.d)
hca.XLabel.String = 'E (eV)';
hca.YLabel.String = sprintf('Differential Particle Flux (%s)',hplusOmni1.units);
legend(hca,'FPI','HPCA','EIS','FEEPS','location','best')
hca.XGrid = 'on';
hca.YGrid = 'on';

hca = h2(2);
loglog(hca,omni_fpi.e,omni_fpi.e.*omni_fpi.d,omni_hpca.e,omni_hpca.e'.*omni_hpca.d,omni_eis.e,omni_eis.e.*omni_eis.d,omni_feeps.e,omni_feeps.e'.*omni_feeps.d)
hca.XLabel.String = 'E (eV)';
hca.YLabel.String = sprintf('Differential Energy Flux (%s)','eV/(cm^2 s sr eV)');
legend(hca,'FPI','HPCA','EIS','FEEPS','location','best')
hca.XGrid = 'on';
hca.YGrid = 'on';