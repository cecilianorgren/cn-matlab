%% Load data
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');
units = irf_units;
ic = 1;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); % torbert, psbl + EDR + aftermath
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:34:30.00Z'); % psbl + EDR
tint_ti = irf.tint('2017-07-11T22:33:10.00Z/2017-07-11T22:33:15.00Z'); % where we plot the spectra as a function of energy
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
% c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
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

c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseEperpVixB? = gseE?perp.resample(gseVixB?.time)+gseVixB?; gseEVixB?.name = ''Eperp+VixB'';',ic)
c_eval('gseEperpVexB? = gseE?perp.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''Eperp+VexB'';',ic)

%% Figure
npanels = 7;
nrows = 2;
ncols = 2;
%h = irf_plot(npanels);
[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.5,'vertical');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z},'comp');
  hca.YLabel.String = {'B (GSE)','(nT)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 0 % gseE
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.x,gseE1.y,gseE1.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % gseE, lowpass
  hca = irf_panel('E lowpass');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{plotE.x,plotE.y,plotE.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 1 % gseVixB
  hca = irf_panel('vixB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVixB1.x,-gseVixB1.y,-gseVixB1.z},'comp');
  hca.YLabel.String = {'-v_i\times B (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % gseEperp+vixB
  hca = irf_panel('E+vixB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseEperpVixB1.x,gseEperpVixB1.y,gseEperpVixB1.z},'comp');
  hca.YLabel.String = {'E_\perp+v_i\times B','(mV/m)'};
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
if 1 % iPDist deflux omni  
  hca = irf_panel('i DEF omni');  
  irf_spectrogram(hca,iPDist1.tlim(tint).deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 0 % Ti part
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
if 0 % iPDist dpflux omni feeps 
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
if 0 % iPDist dpflux omni EIS
  hca = irf_panel('i PEF EIS omni');  
  irf_spectrogram(hca,eisi1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};  
end
if 0 % iPDist dpflux omni  FPI
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
if 0 % iPDist dpflux omni, HPCA  
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
if 1 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  irf_spectrogram(hca,ePDist1.tlim(tint).deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTe1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 1 % Ti, Te
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

irf_plot_axis_align
irf_pl_mark(h,tint_ti)
irf_zoom(h,'x',tint)

% ion
omni_fpi.e = iPDist1.tlim(tint_ti).depend{1}(1,:);
omni_fpi.d = nanmean(iPDist1.tlim(tint_ti).dpflux.omni.data,1);
omni_hpca.e = hplusOmni1.depend{1}(:);
omni_hpca.d = nanmean(hplusOmni1.tlim(tint_ti).data,1);
%omni_eis.e = Omnifluxion_epd_eis_brst_l2.depend{1}(1,:);
%omni_eis.d = nanmean(Omnifluxion_epd_eis_brst_l2.data,1)*1e-3; % 1/keV -> /eV
omni_feeps.e = feepsi1.tlim(tint_ti).depend{1}(:);
omni_feeps.d = nanmean(feepsi1.tlim(tint_ti).data,1)*1e-3; % 1/keV -> /eV

ti_av = mean(gseTi1.tlim(tint_ti).trace.data/3,1);
vi_av = mean(gseVi1.tlim(tint_ti).abs.data,1);
Evi_av = units.mp*(vi_av*1e3).^2/2/units.eV;

% electrons
omnie_fpi.e = ePDist1.tlim(tint_ti).depend{1}(1,:);
omnie_fpi.d = nanmean(ePDist1.tlim(tint_ti).dpflux.omni.data,1);
%omni_eis.e = Omnifluxion_epd_eis_brst_l2.depend{1}(1,:);
%omni_eis.d = nanmean(Omnifluxion_epd_eis_brst_l2.data,1)*1e-3; % 1/keV -> /eV
%omni_feeps.e = feepsi1.tlim(tint_ti).depend{1}(:);
%omni_feeps.d = nanmean(feepsi1.tlim(tint_ti).data,1)*1e-3; % 1/keV -> /eV

te_av = mean(gseTe1.tlim(tint_ti).trace.data/3,1);
ve_av = mean(gseVe1.tlim(tint_ti).abs.data,1);
Eve_av = units.me*(ve_av*1e3).^2/2/units.eV;


%% Panels
isub = 1;
if 1 % DPF ion
  hca = h2(isub); isub = isub + 1;
  loglog(hca,omni_fpi.e,omni_fpi.d,omni_hpca.e,omni_hpca.d,omni_feeps.e,omni_feeps.d) %,omni_eis.e,omni_eis.d
  hca.XLabel.String = 'E (eV)';
  hca.YLabel.String = {'Differential Particle Flux',sprintf('(%s)',hplusOmni1.units)};
  legend(hca,'FPI','HPCA','FEEPS','location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim(1) = 10;
end
if 1 % DEF ion
  hca = h2(isub); isub = isub + 1;
  loglog(hca,omni_fpi.e,omni_fpi.e.*omni_fpi.d,omni_hpca.e,omni_hpca.e'.*omni_hpca.d,omni_feeps.e,omni_feeps.e'.*omni_feeps.d) %,omni_eis.e,omni_eis.e.*omni_eis
  hold(hca,'on')
  loglog(hca,ti_av*[1 1],hca.YLim,'--k')
  loglog(hca,Evi_av*[1 1],hca.YLim,'--','color',[0.8 0.8 0.8])
  hold(hca,'off')
  hca.XLabel.String = 'E (eV)';
  hca.YLabel.String = {'Differential Energy Flux','(eV/(cm^2 s sr eV)'};
  legend(hca,'FPI','HPCA','FEEPS','location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim(1) = 10;
end
if 1 % DEF electron
  hca = h2(isub); isub = isub + 1;
  loglog(hca,omnie_fpi.e,omnie_fpi.e.*omnie_fpi.d)
  hold(hca,'on')
  loglog(hca,te_av*[1 1],hca.YLim,'--k')
  loglog(hca,Eve_av*[1 1],hca.YLim,'--','color',[0.8 0.8 0.8])
  hold(hca,'off')
  hca.XLabel.String = 'E_e (eV)';
  hca.YLabel.String = {'Differential Energy Flux','(eV/(cm^2 s sr eV)'};
  legend(hca,{'FPI electron'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim(1) = 40;
  hca.XTick = 10.^(0:6);
end
%for ip = 1:(nrows*ncols)
%  h2(ip).XLim(1) = 10;
%end