% ions_from_north
% 
units = irf_units;
c_eval('iDEFomniV? = iDEFomni64_?; iDEFomniV?.f = sqrt(iDEFomniV?.f*units.eV*2/units.mp)/1000;')
c_eval('eDEFomniV? = eDEFomni64_?; eDEFomniV?.f = sqrt(eDEFomniV?.f*units.eV*2/units.me)/1000;')
%%
tint = irf.tint('2015-10-16T10:32:50.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere
h = irf_plot(5);
ic = 1;
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.98 0.9],'fontsize',12);
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hout,irf_colormap('space'))  
end
if 1 % i DEF omni 64 
  hca = irf_panel(irf_ssub('i DEF omni 64 ?',ic));  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hout,irf_colormap('space'))   
  try
    hold(hca,'on')
    irf_plot(hca,iEnergyLevel,'k')
    %irf_plot(hca,iEnergyLevelResamp,'k')
    hold(hca,'off')
  end
end
if 1 % i DEF omni 64 MMS3  
  hca = irf_panel(irf_ssub('i DEF omni 64 V ?',ic));  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomniV?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};  
  set(hca,'yscale','lin');
  set(hca,'ytick',[0:500:3000]);  
  colormap(hout,irf_colormap('space'))  
  try
    hold(hca,'on')
    irf_plot(hca,iVelocityLevel,'k')
    %irf_plot(hca,iEnergyLevelResamp,'k')
    hold(hca,'off')
  end
end

irf_zoom(h,'x',tint)
%irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%%
tE = ginput(2);
t1 = tE(1,1);
t2 = tE(2,1);
E1 = tE(1,2);
E2 = tE(2,2);
hcf = gcf;
t0 = hcf.UserData.t_start_epoch;

v1 = cn_eV2v(E1,'eV')/sqrt(1836);
v2 = cn_eV2v(E2,'eV')/sqrt(1836);
x = v1*v2/(v1-v2)*(t2-t1);
x = E1*E2/(E1-E2)*(t2-t1);
%
tt = irf_time([t1;t2]+t0,'epoch>epochtt');
iEnergyLevel = irf.ts_scalar(tt,[E1;E2]);
iVelocityLevel = irf.ts_scalar(tt,[E1;E2]);
%iEnergyLevelResamp = iEnergyLevel.resample(disDist1.time,'linear');