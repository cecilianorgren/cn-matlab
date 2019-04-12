%% Wave speeds, trapping, range 
% burst interval
tint = irf.tint('2017-07-06T08:16:03.00Z/2017-07-06T08:18:13.00Z');
% shorter time interval for dispersion analysis
Tints = irf.tint('2017-07-06T08:16:38.40Z/2017-07-06T08:16:39.00Z'); % second "good" batch
%Tints = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:17.00Z'); % first "good" batch

ic = 1:4;
units = irf_units;
localuser = datastore('local','user');

%% Load data

c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('matB? = gseB?.abs.data;',ic)
if 0
c_eval('matParTe? = facTe?.xx.resample(gseB?.time).data;',ic)
c_eval('matParTi? = facTi?.xx.resample(gseB?.time).data;',ic)
c_eval('matPerTe? = (facTe?.yy.resample(gseB?.time).data + facTe?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matPerTi? = (facTi?.yy.resample(gseB?.time).data + facTi?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matTe? = facTe?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matTi? = facTi?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matNe? = ne?.resample(gseB?.time).data;',ic)
c_eval('vte?perp = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Vte''); vte?perp = irf.ts_scalar(gseB?.time,vte?perp)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte?par = irf_plasma_calc(matB?,matNe?,0,matParTe?,matTi?,''Vte''); vte?par = irf.ts_scalar(gseB?.time,vte?par)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Vte''); vte? = irf.ts_scalar(gseB?.time,vte?)*1e-3; vte?.units = ''km/s'';',ic)
end

disp('Done loading data.')

fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
dirNameMatlab = sprintf('+mms_%s%s%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
eventPath = ['/Users/' localuser '/Research/Events/' dirName '/'];
matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/' dirNameMatlab '/'];

%% % Load wave phase velocities (the ones obtained semi-manually)
fid = fopen([matlabPath 'esw_properties_redo.txt'],'r');
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

if 0 % load esw_properties_irf_4_v_gui.txt
  fid = fopen([matlabPath 'esw_properties_irf_4_v_gui.txt'],'r');
  data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
  esw_data_ = textscan(fid,[data_format_read]);
  fclose(fid)
  for icell = 1:numel(esw_data_)
    esw_data{icell} = cat(1,esw_data{icell},esw_data_{icell});
  end
end

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);    
  esw_data{icell}(8) = [];
end

tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
c_eval('tsVtrap? = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{9+?}/units.me)*1e-3);',1:4)
c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)
vmax = 0.25*(vmax1 + vmax2 + vmax3 + vmax4); vmax.name = 'av(vph+vtrap)';
vmin = 0.25*(vmin1 + vmin2 + vmin3 + vmin4); vmin.name = 'av(vph-vtrap)';
  
%% Dispersion analysis
tic
[xvecs,yvecs,Power] = mms.fk_powerspec4SC('gseE?par','gseR?','gseB?',Tints,'linear',10,'numk',500,'cav',4,'wwidth',2);
toc
% 202 s

%% Polarization analysis, not necessary
tint_polan = Tints + [-1.1 1.1];
ic = 1;
c_eval('Rxyz = gseR?;',ic);
c_eval('Bxyz = gseB?.tlim(tint_polan);',ic);
c_eval('Bscm = gseB?scm.tlim(tint_polan);',ic);
c_eval('Exyz = gseE?.tlim(tint_polan);',ic);


Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
B_SI=Bxyz.abs.data*1e-9;
Wce = e*B_SI/Me;
ecfreq = Wce/2/pi;
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = irf.ts_scalar(Bxyz.time,ecfreq);
ecfreq01 = irf.ts_scalar(Bxyz.time,ecfreq01);
ecfreq05 = irf.ts_scalar(Bxyz.time,ecfreq05);

polarization = irf_ebsp(Exyz,Bscm,Bxyz,Bxyz,Rxyz,[10 4000],'polarization','fac');
%polarization = irf_ebsp(gseE1,gseB1scm,gseB1,gseB1,gseR1,[2 4000],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Bperp = polarization.bb_xxyyzzss(:,:,1)+polarization.bb_xxyyzzss(:,:,2);
Esum = polarization.ee_xxyyzzss(:,:,4);
Eperp = polarization.ee_xxyyzzss(:,:,1)+polarization.ee_xxyyzzss(:,:,2);
Epar = polarization.ee_xxyyzzss(:,:,3);
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);
% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;
vphperp = sqrt(Eperp./Bperp)*1e6;

% Remove points with very low B amplitutes
Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
pfluxz(removepts) = NaN;
vph(removepts) = NaN;
vphperp(removepts) = NaN;

%% Reduced electron distribution
eint = [00 40000];
lowerelim = 100;
vint = [-Inf Inf];
tint_fred = Tints + [-1.1 1.1];
c_eval('eDist = ePDist?.tlim(tint_fred);',ic)
c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
nMC = 500;
tic; ef1D = eDist.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC); toc % reduced distribution along B

%% Plot timeseries and 4SC disprel
npanels = 3;
[h,h2] = initialize_combined_plot(npanels,2,1,0.5,'vertical'); % horizontal

for iip2 = 1:numel(h2)
  h2(iip2).Position(1) = h2(iip2).Position(1) + 0.05;
end

ic = 1;
isub = 0;
zoomy = [];
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % Epar spectra
  isub = isub + 1;
  hca=irf_panel('Eperp spectra');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=0.5*Eperp;
  specrec.f_label='';
  specrec.p_label={'log_{10}E_{\perp}^{2}','mV^2 m^{-2} Hz^{-1}'};
  [~,hcb] = irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
end
if 0 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{tsPhi},'.');
  hca.YLabel.String = {'\phi_{||}','(V)'};  
  hca.YLabel.Interpreter = 'tex';
end

if 0 % v phase
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{tsVphpar},'comp');
  hca.YLabel.String = {'v_{ph}','(km/s)'};  
end
if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');    
  hold(hca,'on')
  vscale = 1e-3;
  
  set(hca,'ColorOrder',mms_colors('111223344'))
  irf_plot(hca,{tsVphpar*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  hca.YLim = sort([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]);
  hold(hca,'off')
  hca.YLabel.String = {'v','(10^3 km/s)'};  
  set(hca,'ColorOrder',mms_colors('122'))
  irf_legend(hca,{'v_{ph}','v_{trap}'},[0.98 0.9],'fontsize',12);    
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('efred');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')

  %set(hca,'ColorOrder',mms_colors('122'))
  %irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  %hca.YLabel.String = 'v_e (km/s)'; 
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
%  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('efred + vtrap');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  vscale = 1e-3;
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))
  irf_plot(hca,{tsVphpar.tlim(Tints)*vscale,vmin.tlim(Tints)*vscale,vmax.tlim(Tints)*vscale},'comp');
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))
  plot(hca,hca.XLim,[0 0],'--','linewidth',1);
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  %hca.YLabel.String = 'v_e (km/s)'; 
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(lowerelim,'%g') ' eV'],[0.01 0.99],'color',0*[1 1 1])
  hca.YLim = vscale*sort([max(real(vmax.tlim(Tints).data)) min(real(vmin.tlim(Tints).data))]);
  
  %hca.CLim = [-5 -2.5];
  %hca.CLim = log10(prctile(specrec.p(:),[10 90]));
end

%irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));
irf_zoom(h,'x',Tints + [-1 1]);
irf_zoom(h(zoomy),'y')
%irf_plot_axis_align(h)
%
% dispersion relation
isub = 1;
if 1
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  vph = mms.estimate_phase_speed(Power.Powerkmagf,yvecs.fkmag,xvecs.kmag,100);
  kfit = 2*pi*yvecs.fkmag/vph;
  hold(hca,'on')
  plot(hca,kfit,yvecs.fkmag*1e-3,'k','linewidth',1)
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 1 % f vx |k|, with all the semi-manual velocities added.
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  
  vph_tint = tsVph.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs.kmag([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
  end
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kperp,yvecs.kpar,log10(Power.Powerkperpkpar));
  shading(hca,'flat');
  xlabel(hca,'|k_{\perp}| (m^{-1})');
  ylabel(hca,'k_{||} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_{\perp},k_{||})/P_{max}');
  colormap('jet');
  irf_legend(hca,'(b)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxkxky,yvecs.kykxky,log10(Power.Powerkxky));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kykxky)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kxkxky)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'k_{y} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,k_y)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(c)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxkxkz,yvecs.kzkxkz,log10(Power.Powerkxkz));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kzkxkz)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kxkxkz)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'k_{z} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,k_z)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(d)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kykykz,yvecs.kzkykz,log10(Power.Powerkykz));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kzkykz)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kykykz)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{y} (m^{-1})');
  ylabel(hca,'k_{z} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_y,k_z)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(e)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxf,yvecs.fkxf*1e-3,log10(Power.Powerkxf));
  shading(hca,'flat');
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(f)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kyf,yvecs.fkyf*1e-3,log10(Power.Powerkyf));
  shading(hca,'flat');
  xlabel(hca,'k_{y} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_y,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(g)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kzf,yvecs.fkzf*1e-3,log10(Power.Powerkzf));
  shading(hca,'flat');
  xlabel(hca,'k_{z} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_z,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(h)',[0.99 0.99],'color','k','fontsize',14)
end

%% Plot timeseries and 4SC disprel
npanels = 1;
nrows = 1;
ncols = 3;
%h = subplot(nrows,ncols,[1 2]);
%h2(1) = subplot(nrows,ncols,3);

h = subplot('position',[0.07 0.15 0.6 0.8]);
h2 = subplot('position',[0.74 0.15 0.23 0.8]);
%h2(2) = subplot(nrows,ncols,3);

fontsize = 12;
ic = 1;
isub = 1;
zoomy = [];

if 1 % e psd vpar  
  hca = h(isub); isub = isub + 1;
  position = hca.Position;
  vscale = 1e-3;
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;  
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))
  if 0
    irf_plot(hca,{tsVphpar.tlim(Tints)*vscale,vmin.tlim(Tints)*vscale,vmax.tlim(Tints)*vscale},'comp');
  else
    linewidth = 2;
    linemarker = '*';
    linestyle = 'none';
    irf_plot(hca,tsVphpar.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker,'linestyle',linestyle);
    irf_plot(hca,vmin.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
    irf_plot(hca,vmax.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
  end
  
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  if 1 % plot v = 0
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('111'))
    plot(hca,hca.XLim,[0 0],'--','linewidth',1);
    hold(hca,'off')
  end
  if 1 % plot vepar
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale,'-');
    hold(hca,'off')
  end
  if 0 % plot vepar +- vth
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+vte1par.resample(gseVe1)*vscale,'-');
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+-1*vte1par.resample(gseVe1)*vscale,'-');
    hold(hca,'off')
  end
  
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';  
  hcb.YLabel.String = 'PSD (s/m^4)';
    
  %hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  %hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(lowerelim,'%g') ' eV'],[0.98 0.05],'color',0*[1 1 1],'fontsize',fontsize)
  hca.YLim = 1*vscale*sort([max(real(vmax.tlim(Tints).data)) min(real(vmin.tlim(Tints).data))]);
  hca.YLim(1) = -30;
  hca.YLim(2) = 20;
  
  hca.CLim = [-5.5 -2.3];
  %hcb.Position(1) = hcb.Position(1) + (position(3) - hca.Position(3));
  %hca.Position = position;
  
end

%irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));
irf_zoom(h,'x',Tints + 1*[-1 1]);

%
% dispersion relation
isub = 1;
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  vph = mms.estimate_phase_speed(Power.Powerkmagf,yvecs.fkmag,xvecs.kmag,100);
  kfit = 2*pi*yvecs.fkmag/vph;
  hold(hca,'on')
  plot(hca,kfit,yvecs.fkmag*1e-3,'k','linewidth',1)
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 1 % f vx |k|, with all the semi-manual velocities added.
  hca = h2(isub); isub = isub + 1;
  kscale = 1e3;
  pcolor(hca,xvecs.kmag*kscale,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'|k| (km^{-1})');
  else
    xlabel(hca,'|k| (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs.kmag([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec*kscale,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
  end
  hold(hca,'off')
  axis(hca,[0 2.3e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
end

h.FontSize = fontsize;
h2.FontSize = fontsize;
