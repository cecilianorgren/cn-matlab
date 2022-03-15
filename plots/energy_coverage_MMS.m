%% Load data
%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 2;

%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% Time from time interval
%tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');
tint_action = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');

% Time from file name
fileId = '20170725220853';

iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));
tint = [files(iFile-1).start files(iFile).stop] + [1 -1];

% Event path
eventPath = ['/Users/' localuser '/Research/Events/mms_' fileId '/']; % for saving figures
%matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms_' fileID '/'];
mkdir(eventPath)

%%
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
% Remove all one-count "noise"
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_remnoise = iPDist?;',ic)
c_eval('iPDist?_remnoise.data(iPDist?_remnoise.data < iPDistErr?.data*1.01) = 0;',ic)

c_eval('eis_omni? = mms.get_data(''Omnifluxproton_epd_eis_brst_l2'',tint,?);',ic)
c_eval('feeps_ion_omni? = mms.get_data(''Omnifluxion_epd_feeps_brst_l2'',tint,?);',ic)

%% Figure
ic = 2;
c_eval('iPDist = iPDist?_remnoise;',ic)

npanels = 3;
%h = irf_plot(npanels);
[h1,h2] = initialize_combined_plot(3,1,1,0.65,'vertical');
iisub = 0;
c_eval('h1(?).Position(2) = h1(?).Position(2) + 0.03;',1:3)
c_eval('h1(?).Position(1) = 0.10;',1:3)
h2.Position(1) = 0.7;
axis(h2,'square')


cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 1 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('specrec = eis_omni?.specrec(''energy'');',ic)
  % add a nan energy level to get the yscale in 10^x format
  %specrec.p = [specrec.p, nan(size(specrec.t))];
  %specrec.f = [specrec.f, 10e4*ones(size(specrec.t))];
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  %hca.YLim(2) = specrec.f(1,end-1);
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',ic)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex'; 
end
if 1 % i DPF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  specrec = iPDist.omni.dpflux.specrec;
  specrec.p = specrec.p*1000; % 1/eV - > 1/keV
  specrec.p_label = {'1/(cm^2 s sr keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log'); 
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i^{FPI}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end

irf_zoom(h1,'x',tint)
irf_plot_axis_align(h1)

if 1 % all instruments averaged and combined
  %%
  isub = 1;
  hca = h2(isub); iub = isub + 1;
  fpi.energy = iPDist.depend{1}(1,:);
  fpi.data = mean(iPDist.omni.dpflux.data,1)*1000;  % 1/eV -> 1/keV  
  fpi.noise = mean(iPDistErr2.omni.dpflux.data,1)*1000;  % 1/eV -> 1/keV  
  eis.energy = eis_omni2.depend{1}(1,:);
  eis.data = mean(eis_omni2.data,1);
  feeps.energy = feeps_ion_omni2.depend{1}(1,:);
  feeps.data = mean(feeps_ion_omni2.data,1);
  
  hl = loglog(hca,fpi.energy,fpi.data,...
                  eis.energy,eis.data,...
                  feeps.energy,feeps.data);
  c_eval('hl(?).LineWidth = 1.5;',1:3)
  if 0
    hold(hca,'on')
    hlerr = loglog(hca,fpi.energy,fpi.noise,'k--','linewidth',1);  
    hold(hca,'off')
  end
  legend(hca,{'FPI','EPD-EIS','EPD-FEEPS'},'box','off')
  hca.YLabel.String = 'Differential particle flux (1/(cm^2 s sr keV))';
  hca.XLabel.String = 'E_i (eV)';
  axis(hca,'square')
  hca.XMinorGrid = 'off';
  hca.XGrid = 'on';
  hca.YMinorGrid = 'off';
  hca.YGrid = 'on';
  hca.Title.String = 'Time-averaged';
end

c_eval('h1(?).FontSize = 12;',1:3)
c_eval('h2(?).FontSize = 12;',1)
%annotation('textarrow',[0.685 0.685]-0.025,[0.47 0.495]+0.00,'string',{'EDI range'},'fontsize',12,'horizontalalignment','center');

%% Annotationsº
%annotation('textarrow',[0.4 0.49],[0.90 .81],'string',{'Energetic ion in reconnection jet'},'fontsize',12,'horizontalalignment','center');
