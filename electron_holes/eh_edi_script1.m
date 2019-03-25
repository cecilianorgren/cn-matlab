
%% Set time interval, database, and save path
ic = 1;
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');

tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z'); % eh edi paper

volume = 'Fountain';
mms.db_init('local_file_db',['/Volumes/' volume '/Data/MMS/']);
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

save_dir_root = ['/Volumes/' volume '/Figures/MMS/EH_EDI/'];
if not(exist(save_dir_root)); mkdir(eventPath); end

%% Load data
tic;
disp('Loading data.');

% Magnetic field
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

% Spacecraft potential
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

% EDI flux
c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)

% Skymap distributions
%c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Remove all one-count "noise"
%c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
%c_eval('ePDist?.data(ePDist?.data<ePDistErr?.data*1.1) = 0;',ic)

% Particle moments
% Pressure
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);

% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

% Velocity
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);

disp('Done loading data.'); 
toc

%% Prepare data
% Parallel electric field
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('facE? = irf_convert_fac(gseE?,gseB?,[1 0 0]);',ic)

% Pseudo potential, int(E)dt, and high-pass filtered
flow = 5;
fhigh = 0;
c_eval('intEdt? = irf_integrate(gseE?par);',ic)
c_eval('intEdt?filt = intEdt?.filt(flow,0,[],5);',ic)

% Power spectrogram of Epar
%c_eval('wavE?par = irf_wavelet(gseE?par.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
%c_eval('wavE?par.f_units = ''Hz''; wavE?fac.f_label = ''f [Hz]''; wavE?fac.p_label = {''log_{10} E_{||}^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavE?fac = irf_wavelet(facE?.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE?fac.f_units = ''Hz''; wavE?fac.f_label = ''f [Hz]''; wavE?fac.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)

% Spectrogram of EDI flux
c_eval('wavFluxEDI? = irf_wavelet(ePitch?_flux_edi.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 500],''nf'',50);',ic)

% Set up specrecs
c_eval('srE?par.p = {wavE?fac.p{3}}; srE?par.t = wavE?fac.t; srE?par.f = wavE?fac.f;',ic)
c_eval('srE?par.f_units = ''Hz''; srE?par.f_label = ''f (Hz)''; srE?par.p_label = {''log_{10} E_{||}^2'',''(mV/m)^2/Hz''};',ic)
c_eval('srE?perp.p = {}; srE?perp.p = {0.5*(wavE?fac.p{1}+wavE?fac.p{2})}; srE?perp.t = wavE?fac.t; srE?perp.f = wavE?fac.f;',ic)
c_eval('srE?perp.f_units = ''Hz''; srE?perp.f_label = ''f (Hz)''; srE?perp.p_label = {''log_{10} E_{\perp}^2'',''(mV/m)^2/Hz''};',ic)
c_eval('srFlux? = wavFluxEDI?;',ic)
c_eval('srFlux?.f_units = ''Hz''; srFlux?.f_label = ''f (Hz)''; srFlux?.p_label = {''log_{10} j^2'',''(cm^{-2}s^{-1}sr^{-1})^2/Hz''};',ic)

% Get CLim for spectrograms
edges_log = -20:1:20;
pdf_lower_lim = 0.01;
percentiles = [1 99];

% EDI Flux
c_eval('climFlux?{!} = prctile(log10(srFlux?.p{!}(:)),percentiles);',ic,1:8)
climFlux = [Inf -Inf];
for iic = 1:1
  for inode = 1:8
    c_eval('climFlux(1) = min([climFlux(1) climFlux?{inode}(1)]);',iic)
    c_eval('climFlux(2) = max([climFlux(2) climFlux?{inode}(2)]);',iic)
  end  
end

% Epar, Eperp
c_eval('climEpar? = prctile(log10(srE?par.p{1}(:)),percentiles);',ic)
c_eval('climEperp? = prctile(log10(srE?perp.p{1}(:)),percentiles);',ic)
climE = [Inf -Inf];
for iic = 1:1
  c_eval('climE(1) = min([climE(1) climEpar?(1) climEperp?(1)]);',iic)
  c_eval('climE(2) = max([climE(2)  climEpar?(2) climEperp?(2)]);',iic)  
end

%% Cross-correlation of intEpar and flux edi
% https://en.wikipedia.org/wiki/Degree_of_coherence#Degree_of_first-order_coherence
tint_test = irf.tint('2017-07-06T13:54:00.00Z/2017-07-06T13:54:10.00Z'); % eh edi paper

c_eval('time_edi? = ePitch?_flux_edi.tlim(tint_test).time;',ic)
c_eval('tE = facE?.resample(time_edi?).z;',ic)
c_eval('tF = ePitch?_flux_edi.tlim(tint_test).resample(time_edi?);',ic)
c_eval('tF = facE?.resample(time_edi!).z;',2,ic)

palim = [168 180];
wE_ = irf_wavelet(tE,'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);
%wF_ = irf_wavelet(tF.palim(palim),'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);
wF_ = irf_wavelet(tF,'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);

wE = wE_.p{1};
wF = wF_.p{1};

pE = wE.*conj(wE);
pF = wF.*conj(wF);

wC = wE.*conj(wF)./sqrt(pE.*pF);
%wC = wE.*conj(wF)./(abs(pE).*abs(pF));

srE.p = abs(wE);      srE.t = wE_.t;     srE.f = wE_.f;     srE.p_label = 'E (abs) (mV/m)';
srEreal.p = real(wE); srEreal.t = wE_.t; srEreal.f = wE_.f; srEreal.p_label = 'E (imag) (mV/m)';
srEimag.p = imag(wE); srEimag.t = wE_.t; srEimag.f = wE_.f; srEimag.p_label = 'E (real) (mV/m)';

srF.p = abs(wF);      srF.t = wF_.t;     srF.f = wF_.f;     srF.p_label = 'F (abs) (...)';
srFreal.p = real(wF); srFreal.t = wF_.t; srFreal.f = wF_.f; srFreal.p_label = 'F (imag) (...)';
srFimag.p = imag(wF); srFimag.t = wF_.t; srFimag.f = wF_.f; srFimag.p_label = 'F (real) (...)';

srC.p = abs(wC);      srC.t = wE_.t;     srC.f = wE_.f;     srC.p_label = 'C (abs)';
srCreal.p = real(wC); srCreal.t = wE_.t; srCreal.f = wE_.f; srCreal.p_label = 'C (imag)';
srCimag.p = imag(wC); srCimag.t = wE_.t; srCimag.f = wE_.f; srCimag.p_label = 'C (real)';


%% Figure: EDI Flux coherence
npanels = 9;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;

step_efw = 50;25;
step_edi = 1;3;

varstrs = {'srE','srEreal','srEimag','srF','srFreal','srFimag','srC','srCreal','srCimag'};
nvars = numel(varstrs);
for ivar = 1:nvars  
  varstr = varstrs{ivar};
  specrec = eval(varstr);  
  specrec.p = specrec.p(1:step_edi:end,:); specrec.t = specrec.t(1:step_edi:end);
  hca = irf_panel(varstr);    
  irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');  
  hca.YTick = 10.^[-10:1:10];
  hca.CLim = sort(real(prctile(log10(specrec.p(:)),[1 99])));
  drawnow;
end
hlink = linkprop(h,{'XLim'});

irf_zoom(h,'x',tint_test)
irf_plot_axis_align
  
%% Figure: Overview
npanels = 5;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');
zoomy = [];

if 1 % e DEF omni
  iisub = iisub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % Te par perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % ne
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % J
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % ePDist pa 64
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels
  h(ii).FontSize = 12;
end

%% Figure: EDI flux, all nodes
npanels = 6;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');
zoomy = [];

step_efw = 50;25;
step_edi = 6;3;

if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 1 % E perp spectrogram
  iisub = iisub + 1;
  hca = irf_panel(sprintf('E perp spectrogram'));  
  c_eval('specrec = srE?perp; specrec.p = specrec.p{1}(1:step_efw:end,:); specrec.t = specrec.t(1:step_efw:end);',ic)
  irf_spectrogram(hca,specrec,'log'); % [hout,hcb] = 
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hold(hca,'on')
  %c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  %hold(hca,'off')
  %hca.YLabel.String = {'f_','()'};   
  %colormap(hca,cmap) 
  hca.YTick = 10.^[-10:1:10];
  hca.CLim = climE;
  drawnow;
end
if 0 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 1 % E par spectrogram
  iisub = iisub + 1;
  hca = irf_panel(sprintf('E par spectrogram'));
  c_eval('specrec = srE?par; specrec.p = specrec.p{1}(1:step_efw:end,:); specrec.t = specrec.t(1:step_efw:end);',ic)
  irf_spectrogram(hca,specrec,'log'); % [hout,hcb] = 
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hold(hca,'on')
  %c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  %hold(hca,'off')
  %hca.YLabel.String = {'f_','()'};   
  %colormap(hca,cmap) 
  hca.YTick = 10.^[-10:1:10];
  hca.CLim = climE;
  drawnow;
end
for inode = [1 8] % edi flux spectrogram
  iisub = iisub + 1;
  hca = irf_panel(sprintf('edi flux pa %g',inode));  
  c_eval('specrec = srFlux?; specrec.p = specrec.p{!}(1:step_edi:end,:); specrec.t = specrec.t(1:step_edi:end); specrec.f = specrec.f;',ic,inode)
  irf_spectrogram(hca,specrec,'log'); % [hout,hcb] = 
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hold(hca,'on')
  %c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  %hold(hca,'off')
  %hca.YLabel.String = {'f_','()'};   
  %colormap(hca,cmap) 
  hca.YTick = 10.^[-10:1:10];
  hca.CLim = climFlux;
  drawnow;
end
if 0 % E j coherence
  %%
  iisub = iisub + 1;  
  hca = irf_panel(sprintf('E j coherence 2'));
  c_eval('specrec = srCoh; specrec.p = specrec.p{1}(1:step_edi:end,:); specrec.t = specrec.t(1:step_edi:end);',ic)
  irf_spectrogram(hca,specrec,'log'); % [hout,hcb] = 
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hold(hca,'on')
  %c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  %hold(hca,'off')
  %hca.YLabel.String = {'f_','()'};   
  %colormap(hca,cmap) 
  hca.YTick = 10.^[-10:1:10];
  hca.CLim = climE;
  drawnow;
end
irf_zoom(h,'x',tint)
irf_plot_axis_align
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

irf_zoom(h(1:npanels),'x',tint)
%irf_zoom(h(zoomy),'y')

h(1).Title.String = irf_ssub('MMS ?',ic);
for ii = 1:npanels
  h(ii).FontSize = 12;
end

