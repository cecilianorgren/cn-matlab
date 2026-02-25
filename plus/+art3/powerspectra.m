% Overview plot of electric and magnetic field spectra from EFW, STA and
% WHI.

tbeam = toepoch([2007 08 31 10 17 38.158; 2007 08 31 10 17 39.675])';
tov  = toepoch([2007 08 31 10 14 00.000; 2007 08 31 10 20 00.000])';

ic = 3;
nPanels = 5;
h = irf_plot(nPanels,'newfigure');
isub=1;

if 1 % B C3 GSM (1 panel)       
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(gsmB3));
    ylabel(hca,'B_{GSM} [nT] ');
    irf_legend(hca,{'x','y','z','|B|'},[0.98 0.92]);
    irf_zoom(hca,'y')
end
if 1 % E C3 ISR2 (1 panel)
    diE3s=diE3;
    diE3s(:,3)=diE3s(:,3)-50;
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3s(:,1:3));
    ylabel(hca,'E_{ISR2} [mV/m]');
    irf_legend(hca,{'x','y'},[0.98 0.92]);
    %irf_legend(h(k),{'C3','C4'},[0.98 0.95]);
end
if 0,   % PANEL: C?       STAFF spectrogram Bx
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('BB_xxyyzz_sr2__C?_CP_STA_PSD',ic);
    irf_plot(hca,varname,'colorbarlabel',c_caa_var_get(varname,'units'),'fitcolorbarlabel','comp',1);
    if 1
    hold(hca,'on');    
    c_eval('fce=irf_plasma_calc(gsmB?,0,0,0,0,''Fce'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
    irf_plot(hca,[fce(:,1) fce(:,2)*.5],'-','linewidth',0.2,'color','w'); % fce/2 white line
    irf_plot(hca,[fce(:,1) fce(:,2)*.25],'-','linewidth',0.2,'color','w'); % fce/4 white line
    c_eval('flh=irf_plasma_calc(gsmB?,1,0,0,0,''Flh'');',ic); % calculate lower hybrid frequency (underdense case in puter magnetosphere)
    irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
    hold(hca,'off');
    end   
    
    caxis(hca,[-10 -7]);
    set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
    irf_zoom(hca,'y',[50 4000]);
end
if 0,   % PANEL: C?       STAFF spectrogram Ex and fce/flh lines
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('EE_xxyy_sr2__C?_CP_STA_PSD',ic);
    %varunits=c_caa_var_get(varname,'unit');
    varunits='(mV/m)^2/Hz';
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','comp',1,'nolabels');
    % next lines are examples how to add gyro/lower hybrid frequency lines
    % B & n should be calculated before
    hold(hca,'on');
    if 1
    c_eval('fce=irf_plasma_calc(gsmB?,0,0,0,0,''Fce'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
    c_eval('fpe=irf_plasma_calc(irf_resamp(gsmB?,peaNe?),peaNe?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fpe,'-','linewidth',0.2,'color','k');
    c_eval('fpemom=irf_plasma_calc(irf_resamp(gsmB?,peaNe?),peaNe?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fpemom,'.','linewidth',0.2,'color','r');
    c_eval('flh=irf_plasma_calc(gsmB?,1,0,0,0,''Flh'');',ic); % calculate lower hybrid frequency (underdense case in puter magnetosphere)
    irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
    end
    hold(hca,'off');
    % polish the panel
    caxis(hca,[-9 -1]);
    set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
    irf_zoom(hca,'y',[50 4000]);
end  
if 0,   % PANEL: C?       WHISPER spectrogram
  hca=h(isub); isub=isub+1;
  varname=irf_ssub('Electric_Spectral_Power_Density__C?_CP_WHI_NATURAL',ic);
  varunits='(V/m)^2/Hz';
  % REMOVE 'fillspectgrogramgaps' flag in the next line if precise intervals of
  % WHISPER measurements are needed !!!!
  % If working with shorter intervals can also remove 'tint' option
  irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','fillspectrogramgaps','nolabels');
  if 1, % add plasma frequency lines
      % density n should be calculated before
      hold(hca,'on');
      c_eval('fpe=irf_plasma_calc(irf_resamp(1,peaNe?),peaNe?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
      irf_plot(hca,irf_tappl(fpe,'/1e3'),'-','linewidth',0.2,'color','k');
      hold(hca,'off');
  end
  % polish the panel
  caxis(hca,[-16 -11]);
  set(hca,'yscale','log','ytick',[3 4 5 1e1 20 30 50 ]);
  irf_zoom(hca,'y',[2 12]);
end

%C3_CP_RAP_ESPCT6
%C3_CP_RAP_I3DM_H
%C3_CP_RAP_PAD_E3DD
%C3_CP_RAP_PEDPOS_BM
%C3_CP_RAP_PED_BM
if 0,   % PANEL: C?       RAPID ion spectrogram
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Proton_Dif_flux__C?_CP_RAP_ESPCT6',ic);
    %varunits=irf_get_data(varname,'caa','units');
    varunits={'log_{10} dF','1/cm^2 s sr keV'};
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[0.51 4.49]);
    set(hca,'yscale','log');
    ylabel(hca,'E [keV]');
    irf_zoom(hca,'y',[80 2000]);
    set(hca,'ytick',[1e2 2e2 5e2 1e3 1e4 1e5])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 1,   % PANEL: C?       RAPID PAD_E3DD spectrogram parallel
	ic = 3;
	hca=h(isub); isub=isub+1;
	varname=irf_ssub('PAD_Electron_Dif_flux__C?_CP_RAP_PAD_E3DD',ic);
	irf_plot(hca,varname,'colorbarlabel',varunits,...
		'fitcolorbarlabel','comp',1);
	caxis(hca,[1.1 4.3]);
	set(hca,'yscale','log');
	set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5]);
end
if 0,   % PANEL: C?       RAPID spectrogram pitch angle
	hca=h(isub); isub=isub+1;
	varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
	varunits=irf_get_data(varname,'caa','units');
	irf_plot(hca,varname,'colorbarlabel',varunits,...
		'fitcolorbarlabel','comp_dim1','comp',1);
	caxis(hca,[1.1 4.3]);
	set(hca,'yscale','lin');
	set(hca,'ytick',[0 45 90 135 180]);
	ylabel(hca,'Pitch ang. [deg]');
end
irf_zoom(h,'x',tov)
irf_pl_mark(h,tbeam)
irf_plot_axis_align;


