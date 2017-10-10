% Plots E and B time series and of burst mode electric field in GSE
% coordinates and field-aligned coordinates. Plots spectrograms of parallel
% and perpendicular electric fields and fluctuating magnetic field. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number

Tint = irf.tint('2015-11-30T00:21:44.00Z/2015-11-30T00:26:44.00Z');
TintWavelet = irf.tint('2015-11-30T00:24:20.00Z/2015-11-30T00:24:30.00Z');
%% Load Data 
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',Tint);',ic);
magB = Bxyz.abs;
Bxyzmag = TSeries(Bxyz.time,[Bxyz.data magB.data]);

%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);

%% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm.time.epochUnix));
Exyzfachf = Exyzfac.filt(10,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
Bscmfachf = Bscmfac.filt(10,0,dfB,5);

%% Wavelet transforms
if 1
  nf = 100;
  tic; Ewavelet = irf_wavelet(Exyzfac.tlim(TintWavelet),'nf',nf,'f',[5 4000]); toc % 3 min = 140 s
  tic; Bwavelet = irf_wavelet(Bscmfac.tlim(TintWavelet),'nf',nf,'f',[5 4000]); toc % 3 min = 140 s
else
  nf = 100;
  load /Users/Cecilia/Data/MMS/20151130_wavelets.mat
  c_eval('Ewavelet = wavE?;',ic)
  c_eval('Bwavelet = wavB?;',ic)
end
%%
%compress wavelet transform data into nc-point average
nc = 10;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specperpE=struct('t',Ewavelettimes);
specperpE.f=Ewavelet.f;
specperpE.p=Ewaveletx+Ewavelety;
specperpE.f_label='';
specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};

specparE=struct('t',Ewavelettimes);
specparE.f=Ewavelet.f;
specparE.p=Ewaveletz;
specparE.f_label='';
specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};


idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
Bwavelettimes = Bwavelet.t(idx);
Bwaveletx = zeros(length(idx),nf);
Bwavelety = zeros(length(idx),nf);
Bwaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
        Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specperpB=struct('t',Bwavelettimes);
specperpB.f=Bwavelet.f;
specperpB.p=Bwaveletx+Bwavelety;
specperpB.f_label='';
specperpB.p_label={'log_{10} B_{\perp}^2','nT^2 Hz^{-1}'};

specparB=struct('t',Bwavelettimes);
specparB.f=Bwavelet.f;
specparB.p=Bwaveletz;
specparB.f_label='';
specparB.p_label={'log_{10} B_{||}^2','nT^2 Hz^{-1}'};

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);

%% Plot Figure

h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

if 0
xwidth = 0.86;
ywidth = 0.13;
set(h(1),'position',[0.10 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.10 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.10 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.10 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.10 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.10 0.97-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.10 0.97-7*ywidth xwidth ywidth]);
end

hca=irf_panel('Bxyz');
irf_plot(hca,Bxyzmag);
ylabel(hca,{'B','(nT)'},'Interpreter','tex');
%irf_zoom(hca,'y',[-50 60]);
irf_legend(hca,{'B_{x}','B_{y}','B_{z}','|B|'},[0.98 0.12])
irf_legend(hca,'(a)',[0.99 0.94],'color','k','fontsize',12)

hca=irf_panel('Elf');
irf_plot(hca,Exyzfaclf);
ylabel(hca,{'E (mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.98 0.12])
irf_legend(hca,'(b)',[0.99 0.94],'color','k','fontsize',12)

hca=irf_panel('Ehf');
irf_plot(hca,Exyzfachf);
ylabel(hca,{'\delta E (mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.98 0.12])
irf_legend(hca,'(c)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(hca,'f > 10 Hz',[0.1 0.1],'color','k','fontsize',12)

hca=irf_panel('Especperp');
irf_spectrogram(h(4),specperpE,'log');
hold(hca,'on');
irf_plot(hca,Flh,'color','k','LineWidth',1.5)
irf_plot(hca,Fce,'color','r','LineWidth',1.5)
irf_plot(hca,Fpp,'color','b','LineWidth',1.5)
hold(hca,'off');
irf_legend(hca,'(d)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(hca,'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(hca,'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(hca,'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(hca,[-6 0]);
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(5)=irf_panel('Especpar');
irf_spectrogram(h(5),specparE,'log');
hold(h(5),'on');
irf_plot(h(5),Flh,'color','k','LineWidth',1.5)
irf_plot(h(5),Fce,'color','r','LineWidth',1.5)
irf_plot(h(5),Fpp,'color','b','LineWidth',1.5)
hold(h(5),'off');
irf_legend(h(5),'(e)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(h(5),'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(h(5),'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(h(5),'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(h(5),[-6 0]);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(6)=irf_panel('Bscmhf');
irf_plot(h(6),Bscmfachf);
ylabel(h(6),{'\delta B (nT)'},'Interpreter','tex');
irf_legend(h(6),{'B_{\perp 1}','B_{\perp 2}','B_{||}'},[0.98 0.12])
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(6),'f > 10 Hz',[0.1 0.1],'color','k','fontsize',12)

hca=irf_panel('Bspec perp');
irf_spectrogram(hca,specperpB,'log');
hold(hca,'on');
irf_plot(hca,Flh,'color','k','LineWidth',1.5)
irf_plot(hca,Fce,'color','r','LineWidth',1.5)
irf_plot(hca,Fpp,'color','b','LineWidth',1.5)
hold(hca,'off');
irf_legend(hca,'(g)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(hca,'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(hca,'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(hca,'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(hca,[-8 -1]);
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'f','(Hz)'},'fontsize',12,'Interpreter','tex');

hca=irf_panel('Bspec par');
irf_spectrogram(hca,specparB,'log');
hold(hca,'on');
irf_plot(hca,Flh,'color','k','LineWidth',1.5)
irf_plot(hca,Fce,'color','r','LineWidth',1.5)
irf_plot(hca,Fpp,'color','b','LineWidth',1.5)
hold(hca,'off');
irf_legend(hca,'(g)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(hca,'f_{LH}',[0.2 0.60],'color','k','fontsize',12)
irf_legend(hca,'f_{ce}',[0.15 0.60],'color','r','fontsize',12)
irf_legend(hca,'f_{pi}',[0.25 0.60],'color','b','fontsize',12)
caxis(hca,[-8 -1]);
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'f','(Hz)'},'fontsize',12,'Interpreter','tex');


load('caa/cmap.mat');
colormap(h(4),cmap);
colormap(h(5),cmap);
colormap(h(7),cmap);
colormap(h(8),cmap);

c_eval('title(h(1),''MMS?'')',ic);

irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
irf_zoom(h([1:3 6]),'y');
set(h(1:7),'fontsize',12);