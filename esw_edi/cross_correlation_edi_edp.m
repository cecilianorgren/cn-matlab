units = irf_units;
irf.log('critical')
ic = 1;

localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   

%%
tint_all = irf.tint('2015-01-01T00:00:00.00Z/2022-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_edi_brst_l2_amb-pm2',tint_all);
nFiles = numel(files);



for iFiles = 1%:nFiles
  %tint = [files(iFile).start files(iFile).stop];
  % EDI
  c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)
  %c_eval('ePitch?_flux_edi_err = mms.get_data(''Flux-err-amb-pm2_edi_brst_l2'',tint,?);',ic)
  %c_eval('ePitch?_flux_edi_err_plus   = ePitch?_flux_edi+(ePitch?_flux_edi_err*+1);',ic)
  %c_eval('ePitch?_flux_edi_err_minus  = ePitch?_flux_edi+(ePitch?_flux_edi_err*-1);',ic)
  %c_eval('ePitch?_counts_edi = mms.get_data(''Counts-amb-pm2_edi_brst_l1a'',tint,?);',ic)
  
  % FGM
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  
  % EDP
  c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
  
  %% Cross correlation
  c_eval('intE? = irf_integrate(gseE?par);',ic)
  c_eval('data1 = intE?.data;',ic)
  c_eval('data2a = ePitch?_flux_edi.palim(180).resample(gseE?).data;',ic)
  c_eval('data2b = ePitch?_flux_edi.palim(0).resample(gseE?).data;',ic)
  Ca = xcorr(data1,data2a,10,'coeff');
  Cb = xcorr(data1,data2b,10,'coeff');
  
  %specrec = irf_powerfft(ePitch1_flux_edi.palim(180), 128, 1000, [20, 16]);
  % Make power spectrograms
  nfft = 512;
  sfreqE = 1/(gseE1.time(2)-gseE1.time(1));
  sfreqFlux = 1/(ePitch1_flux_edi.time(2)-ePitch1_flux_edi.time(1));
  overlap = 0.50;
  c_eval('fftE? = irf_powerfft(intE?,nfft,sfreqE,[overlap]);?',ic)
  c_eval('fftFlux?_0 = irf_powerfft(ePitch?_flux_edi.palim(0),nfft,sfreqFlux,[overlap]);?',ic)
  c_eval('fftFlux?_180 = irf_powerfft(ePitch?_flux_edi.palim(180),nfft,sfreqFlux,[overlap]);?',ic)
  
  
  %% Plot data
  h = irf_plot(6);
  
  if 1 % B gse
    hca = irf_panel('B gse');
    set(hca,'ColorOrder',mms_colors('xyza'))  
    c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
    hca.YLabel.String = {'B_{GSE}','(nT)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  end 
  if 1 % E par
    hca = irf_panel('E par');
    set(hca,'ColorOrder',mms_colors('xyza'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_zoom(hca,'y')
  end
  if 1 % E power spectrogram
    hca = irf_panel('Epar fft');
    c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',ic);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(nT)^2/Hz';    
%     c_eval('irf_plot(hca,flh?)',ic);
%     c_eval('irf_plot(hca,fce?)',ic);
    %hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';
  end
  if 1 % EDI flux par antipar
    hca = irf_panel('Flux EDI');
    set(hca,'ColorOrder',mms_colors('xyza'))
    c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(0),ePitch?_flux_edi.palim(180)},''comp'');',ic)
    hca.YLabel.String = {'j^{EDI}','()'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'0^{\circ}','180^{\circ}'},[0.98 0.9],'fontsize',12);
    irf_zoom(hca,'y')
  end
  if 1 % J flux power spectrogram 0
    hca = irf_panel('Flux EDI fft 0');
    c_eval('specrec = fftFlux?_0;',ic)
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f);
    hca.YScale = 'log';
    hcb = colorbar('peer',hca);
    hcb.Label.String = 'j_{EDI}^{0^{\circ}}';   
    hca.YLabel.String = 'f [Hz]';
  end
  if 1 % J flux power spectrogram 180
    hca = irf_panel('Flux EDI fft 180');
    c_eval('specrec = fftFlux?_180;',ic)
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f);
    hca.YScale = 'log';
    hcb = colorbar('peer',hca);
    hcb.Label.String = 'j_{EDI}^{180^{\circ}}';
    hca.YLabel.String = 'f [Hz]';
  end
  
  legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
  nInd = 1;
  for ii = 1:numel(h)
    irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
    nInd = nInd + 1;
    h(ii).FontSize = 12;
  end
  
  irf_zoom(h,'x',tint)
  irf_zoom(h,'y')
  irf_plot_axis_align
end