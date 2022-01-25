%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 2;

mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% Time from time interval
%tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');

% Time from file name
fileId = '20170725220853';

iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));
tint = [files(iFile-1).start files(iFile).stop] + [1 -1];

% Event path
eventPath = ['/Users/' localuser '/Research/Events/mms_' fileId '/']; % for saving figures
%matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms_' fileID '/'];
mkdir(eventPath)

%% Load and calculate density
% Magnetic field
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
% Spacecraft potential
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
% Distributions_fpi
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
% Distribution errors
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary dat
c_eval('ePDistErr? = mms.get_data(''PDERRe_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
% Copy and remove background (remove all one-count "noise")
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic) % missing some ancillary data?
c_eval('ePDist?_nobg = ePDist?; ePDist?_nobg.data(ePDist?_nobg.data < ePDistErr?.data*1.01) = 0;',ic) % missing some ancillary data?

% Calculate density from scPot
c_eval('[ne?_sc,Iph0,Tph0,Iph1,Tph1] = mms.scpot2ne(scPot?,ne?,gseTe?);',ic)
% Calculate density from FPI VDF
energyrange = [0 40000];
c_eval('mome? = mms.psd_moments(ePDist?,scPot?,''energyrange'',energyrange);',ic)
c_eval('momi? = mms.psd_moments(iPDist?,scPot?,''energyrange'',energyrange);',ic)
% With one-count noise removed
c_eval('mome?_nobg = mms.psd_moments(ePDist?_nobg,scPot?,''energyrange'',energyrange);',ic)
c_eval('momi?_nobg = mms.psd_moments(iPDist?_nobg,scPot?,''energyrange'',energyrange);',ic)
% Original but with lower energy lim
energyrange_lim_e = [70 40000];
energyrange_lim_i = [100 40000];
c_eval('mome?_elim = mms.psd_moments(ePDist?_nobg,scPot?,''energyrange'',energyrange_lim_e);',ic)
c_eval('momi?_elim = mms.psd_moments(iPDist?_nobg,scPot?,''energyrange'',energyrange_lim_i);',ic)

% Calculate wave spectrogram, to see if there is some emission line
tint_wav = irf.tint('2017-07-25T22:07:00.00Z/2017-07-25T22:08:00.00Z');
c_eval('wavE? = irf_wavelet(gseE?.tlim(tint_wav));',ic);
c_eval('samplFreqHz = 1/(gseE?.time(2)-gseE?.time(1));',ic);
nFt = samplFreqHz;
overlapPercent = 0; 
smoothWidth = 0;
c_eval('fftE? = irf_powerfft(gseE?,nFt,samplFreqHz,overlapPercent,smoothWidth);',ic);

% Calculate plasma frequency
c_eval('fpe? = irf.ts_scalar(ne?.time,sqrt(units.e^2*ne?.data*1e6/units.eps0/units.me)/2/pi); fpe?.units = ''Hz'';',ic) % Hz
c_eval('fce? = irf.ts_scalar(gseB?.time,sqrt(units.e*gseB?.abs.data*1e-6/units.me)/2/pi); fce?.units = ''Hz'';',ic) % Hz

%% Figure
ic = 2;

npanels = 11;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 1 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % ne,i FPI
  hca = irf_panel('n fpi');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.98])
  hca.YLabel.String = {'n^{FPI}','(cm^{-3})'};
end
if 1 % ne,i scPot
  hca = irf_panel('n scpot');
  set(hca,'ColorOrder',mms_colors('12'))
  %c_eval('toplot = ne?_sc; toplot.data(toplot.data)',ic)
  c_eval('irf_plot(hca,{ne?_sc},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_{e,scPot}'},[0.98 0.98])
  hca.YLabel.String = {'n_e^{scPot}','(cm^{-3})'};
end

if 1 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % scpot
    hold(hca,'on')    
    c_eval('irf_plot(hca,scPot?,''k'')',ic)
    hold(hca,'off')
  end
  if 1 % elim
    hold(hca,'on')    
    c_eval('irf_plot(hca,irf.ts_scalar(tint,[1 1]''*energyrange_lim_e),''k'')',ic)
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % elim
    hold(hca,'on')    
    c_eval('irf_plot(hca,irf.ts_scalar(tint,[1 1]''*energyrange_lim_i),''k'')',ic)
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % e DEF omni nobg
  isub = isub + 1;
  hca = irf_panel('e DEF omni nobg');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?_nobg.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % scpot
    hold(hca,'on')    
    c_eval('irf_plot(hca,scPot?,''k'')',ic)
    hold(hca,'off')
  end
  if 0 % elim
    hold(hca,'on')    
    c_eval('irf_plot(hca,irf.ts_scalar(tint,[1 1]''*energyrange_lim_e),''k'')',ic)
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,'nobg: one-count noise removed',[0.02 0.1],'k')
end
if 1 % i DEF omni nobg
  isub = isub + 1;
  hca = irf_panel('i DEF omni nobg');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_nobg.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % elim
    hold(hca,'on')    
    c_eval('irf_plot(hca,irf.ts_scalar(tint,[1 1]''*energyrange_lim_i),''k'')',ic)
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';  
  irf_legend(hca,'nobg: one-count noise removed',[0.02 0.1],'k')
end
if 1 % ne,i mms.psd_moments 1
  hca = irf_panel('n psdmoments 1');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mome?.n_psd,momi?.n_psd,},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'n_{e}','n_{i}'},[0.98 0.98])
  hca.YLabel.String = {'n^{mms.psd\_moments}','(cm^{-3})'};
end
if 1 % ne,i mms.psd_moments 2
  hca = irf_panel('n psdmoments 2');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mome?_nobg.n_psd,momi?_nobg.n_psd,mome?_elim.n_psd,momi?_elim.n_psd},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'n_{e}^{nobg}','n_{i}^{nobg}',sprintf('n_{e} (%g<E_e<%g eV)',energyrange_lim_e(1),energyrange_lim_e(2)),sprintf('n_{i} (%g<E_i<%g eV)',energyrange_lim_i(1),energyrange_lim_i(2))}',[1.01 0.98])
  hca.YLabel.String = {'n^{mms.psd\_moments}','(cm^{-3})'};
end

if 1 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E spectrogram
  isub = isub + 1;
  hca = irf_panel('E fft');  
  c_eval('specrec = fftE?;',ic)
  specrec.p{1} = specrec.p{1} + specrec.p{2} + specrec.p{3}; 
  specrec.f_unit = 'Hz';  
  specrec.f_label = 'f (Hz)';  
  specrec.p_label = '((mV)^2/Hz)';  
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %hcb = colorbar('peer',hca);  
  set(hca,'ytick',[1e0 1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % fpe
    %%
    hold(hca,'on')
    c_eval('hfpe = irf_plot(hca,fpe?,''color'',''k'',''linestyle'','':'',''linewidth'',0.5);',ic)
    c_eval('hfce = irf_plot(hca,fce?,''color'',''b'',''linestyle'','':'');',ic)
    hold(hca,'off')
  end
  hca.YLabel.String = {'f','(Hz)'};   
  hca.YLabel.Interpreter = 'tex';   
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end


irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);


