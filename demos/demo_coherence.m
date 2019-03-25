%% Set time interval, database, and save path
volume = 'Fountain';
mms.db_init('local_file_db',['/Volumes/' volume '/Data/MMS/']);
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

save_dir_root = ['/Volumes/' volume '/Figures/MMS/EH_EDI/'];
if not(exist(save_dir_root)); mkdir(eventPath); end

tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z'); % eh edi paper
tint_test = tint;%irf.tint('2017-07-06T13:54:25.00Z/2017-07-06T13:54:50.00Z'); % eh edi paper

%% Load data
ic = 1:2;
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('facE? = irf_convert_fac(gseE?,gseB?,[1 0 0]);',ic)
c_eval('S? = facE?.tlim(tint_test).z;',ic)

palim = [168 180];
%palim = [0 11];

timeline = ePitch1_flux_edi.tlim(tint_test);
S1 = ePitch1_flux_edi.tlim(tint_test).palim(palim).resample(timeline);
S2 = facE1.tlim(tint_test).resample(timeline).z;


fs = 1/(S1.time(2) - S1.time(1));
%wcoherence(S1.data,S2.data,seconds(1/100));


frange = [0 500];
pS1 = irf_wavelet(S1,'wavelet_width',5.36*1 ,'nf',50,'returnpower',1);
pS2 = irf_wavelet(S2,'wavelet_width',5.36*1 ,'nf',50,'returnpower',1);

[WCOH,WCS,F] = wcoherence(S1.data,S2.data,fs);
srC.p = {WCOH'};
srC.t = S1.time.epochUnix;
srC.f = F;

%% Figure: Signal coherence
npanels = 5;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;

step_efw = 1;50;25;
step_edi = 5;

ffilt = [];[20 500];
igroup = 1;
varstrs = {'S1','S2'};
hg = {};
nvars = numel(varstrs);
for ivar = 1:nvars  
  varstr = varstrs{ivar};
  variable = eval(varstr); 
  if not(isempty(ffilt)), variable = variable.filt(ffilt(1),ffilt(2),[],5); end
  hca = irf_panel(varstr);  
  hg{igroup}(ivar) = hca;
  irf_plot(hca,variable);
  hca.YLabel.String = varstr;
end

igroup = 2;
varstrs = {'pS1','pS2','srC'};
pscale = {'log','log','lin'};
nvars = numel(varstrs);
for ivar = 1:nvars  
  varstr = varstrs{ivar};
  specrec = eval(varstr);  
  specrec.p = specrec.p{1}(1:step_edi:end,:); specrec.t = specrec.t(1:step_edi:end);
  hca = irf_panel(varstr);   
  hg{igroup}(ivar) = hca;
  irf_spectrogram(hca,specrec,pscale{ivar});
  set(hca,'yscale','log');  
  hca.YTick = 10.^[-10:1:10];
  %hca.CLim = sort(real(prctile(log10(specrec.p(:)),[1 99])));
  drawnow;
end
hlink = linkprop(h,{'XLim'});
linkprop(hg{2},{'YLim'});

irf_zoom(h,'x',tint_test)
irf_plot_axis_align
  


%% Matlab example
t = 0:0.001:2;
x = cos(2*pi*10*t).*(t>=0.5 & t<1.1)+ ...
    cos(2*pi*50*t).*(t>= 0.2 & t< 1.4)+0.25*randn(size(t));
y = sin(2*pi*10*t).*(t>=0.6 & t<1.2)+...
    sin(2*pi*50*t).*(t>= 0.4 & t<1.6)+ 0.35*randn(size(t));
wcoherence(x,y,1000)

%%

c_eval('time_edi? = ePitch?_flux_edi.tlim(tint_test).time;',ic)
c_eval('tE = facE?.resample(time_edi?).z;',ic)
c_eval('tF = ePitch?_flux_edi.tlim(tint_test).resample(time_edi?);',ic)
c_eval('tF = facE?.resample(time_edi!).z;',2,ic)

palim = [168 180];
wE_ = irf_wavelet(tE,'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);
%wF_ = irf_wavelet(tF.palim(palim),'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);
wF_ = irf_wavelet(tF,'wavelet_width',5.36*2 ,'f',[1 500],'nf',50,'returnpower',0);






