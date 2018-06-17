%%
ic = 4;

interval = '20170706135303';
yyyy = interval(1:4);
mm = interval(5:6);
dd = interval(7:8);
starttime = interval(9:end);

Tint = irf.tint('2017-07-06T13:53:50.00Z/2017-07-06T13:54:20.00Z');


c_eval('tmpDataObj = dataobj(mms.get_filepath(''mms?_fgm_brst_l2'',Tint));',ic)
c_eval('Bxyz = mms.variable2ts(get_variable(tmpDataObj,''mms?_fgm_b_gse_brst_l2''));',ic);
magB = Bxyz.abs;
Bxyzmag = TSeries(Bxyz.time,[Bxyz.data magB.data]);

tmp = ['tmpDataObj = dataobj(mms.get_filepath(''mms?_edp_brst_l2_dce'',Tint));'];
c_eval(tmp,ic);
c_eval('Exyz = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_dce_gse_brst_l2''));',ic);

Efac = irf_convert_fac(Exyz, Bxyz, [1, 0, 0]);
[Epar,Eperp]=irf_dec_parperp(Bxyz,Exyz);

tmp = ['tmpDataObj = dataobj(mms.get_filepath(''mms?_edi_brst_l2_amb-pm2'',Tint));'];
%tmp = ['tmpDataObj = dataobj(''/Volumes/DansHD3/data/mms?/edi/brst/l2/amb-pm2/' yyyy '/' mm '/' dd '/mms?_edi_brst_l2_amb-pm2_' interval '_v3.2.0.cdf'');'];
c_eval(tmp,ic);
c_eval('flux0a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_0_brst_l2''));',ic);
c_eval('flux0b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_0_brst_l2''));',ic);
c_eval('flux180a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_180_brst_l2''));',ic);
c_eval('flux180b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_180_brst_l2''));',ic);

flux0180 = irf.ts_scalar(flux180a.time,[(flux0a.data+flux0b.data)/2 (flux180a.data+flux180b.data)/2]);

c_eval('energy1 = get_variable(tmpDataObj,''mms?_edi_energy_gdu1_brst_l2'');',ic);
c_eval('energy2 = get_variable(tmpDataObj,''mms?_edi_energy_gdu2_brst_l2'');',ic);

EDIen = single(energy1.data(1));

tmp = ['tmpDataObj = dataobj(mms.get_filepath(''mms?_fpi_brst_l2_des-moms'',Tint));'];
%tmp = ['tmpDataObj = dataobj(''/Volumes/DansHD3/data/mms?/fpi/brst/l2/des-moms/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-moms_' interval '_v3.2.0.cdf'');'];
c_eval(tmp,ic);
c_eval('ne = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_numberdensity_brst''));',ic);

Bxyzmag = Bxyzmag.tlim(Tint);
Epar = Epar.tlim(Tint);
magB = magB.tlim(Tint);
flux0180 = flux0180.tlim(Tint);

%%
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(magB).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(magB).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fpe = irf.ts_scalar(magB.time,Fpe);

%%
nf = 100; fmin = 10; fmax = 4000;
nfe = 100; fmine = 10; fmaxe = 500;
Ewavelet = irf_wavelet(Epar,'nf',nf,'f',[fmin fmax]);
fwavelet = irf_wavelet(flux0180,'nf',nfe,'f',[fmine fmaxe]);

nc = 50;
idx = nc/2:nc:length(Ewavelet.t)-nc/2;
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
for ii = 1:length(idx)
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specE=struct('t',Ewavelettimes);
specE.f=Ewavelet.f;
specE.p=Ewaveletx;
specE.f_label='';
specE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};

nc = 10;
idx = nc/2:nc:length(fwavelet.t)-nc/2;
fwavelettimes = fwavelet.t(idx);
fwaveletx = zeros(length(idx),nf);
fwavelety = zeros(length(idx),nf);
for ii = 1:length(idx)
        fwaveletx(ii,:) = squeeze(irf.nanmean(fwavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
        fwavelety(ii,:) = squeeze(irf.nanmean(fwavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specf0=struct('t',fwavelettimes);
specf0.f=fwavelet.f;
specf0.p=fwaveletx;
specf0.f_label='';
specf0.p_label={'log_{10} flux^2 par','cm^4 s^{-1}'};

specf180=struct('t',fwavelettimes);
specf180.f=fwavelet.f;
specf180.p=fwavelety;
specf180.f_label='';
specf180.p_label={'log_{10} flux^2 apar','cm^4 s^{-1}'};

%%
h=irf_plot(6,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.86;
ywidth = 0.145;
set(h(1),'position',[0.08 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.08 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.08 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.08 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.08 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.08 0.97-6*ywidth xwidth ywidth]);


h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyzmag);
ylabel(h(1),{'B (nT)'},'Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}','|B|'},[0.98 0.12])
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',12)

h(2)=irf_panel('Epar');
irf_plot(h(2),Epar);
ylabel(h(2),{'E_{||} (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',12)

h(3)=irf_panel('Especperp');
irf_spectrogram(h(3),specE,'log');
hold(h(3),'on');
irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
irf_plot(h(3),Fpp,'color','b','LineWidth',1.5)
irf_plot(h(3),Fpe,'color','m','LineWidth',1.5)
hold(h(3),'off');
irf_legend(h(3),'(c)',[0.99 0.8],'color','k','fontsize',12)
irf_legend(h(3),'f_{LH}',[0.05 0.60],'color','k','fontsize',12)
irf_legend(h(3),'f_{ce}',[0.1 0.60],'color','r','fontsize',12)
irf_legend(h(3),'f_{pi}',[0.15 0.60],'color','b','fontsize',12)
irf_legend(h(3),'f_{pe}',[0.20 0.60],'color','m','fontsize',12)
caxis(h(3),[-6 2]);
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(3),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(4)=irf_panel('EDIflux');
irf_plot(h(4),flux0180);
ylabel(h(4),{'flux (cm^2 s^{-1})'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(4),{'par','apar'},[0.8 0.94])
irf_legend(h(4),sprintf('E_{EDI} = %.0f eV',EDIen),[0.1 0.9],'color','k','fontsize',12)

h(5)=irf_panel('parfluxspec');
irf_spectrogram(h(5),specf0,'log');
hold(h(5),'on');
irf_plot(h(5),Flh,'color','k','LineWidth',1.5)
irf_plot(h(5),Fce,'color','r','LineWidth',1.5)
irf_plot(h(5),Fpp,'color','b','LineWidth',1.5)
hold(h(5),'off');
irf_legend(h(5),'(e)',[0.99 0.8],'color','k','fontsize',12)
caxis(h(5),[7 11]);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

h(6)=irf_panel('aparfluxspec');
irf_spectrogram(h(6),specf180,'log');
hold(h(6),'on');
irf_plot(h(6),Flh,'color','k','LineWidth',1.5)
irf_plot(h(6),Fce,'color','r','LineWidth',1.5)
irf_plot(h(6),Fpp,'color','b','LineWidth',1.5)
hold(h(6),'off');
irf_legend(h(6),'(f)',[0.99 0.8],'color','k','fontsize',12)
caxis(h(6),[7 11]);
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),{'f (Hz)'},'fontsize',12,'Interpreter','tex');

colormap(h(3),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');

c_eval('title(h(1),''MMS?'')',ic);

irf_plot_axis_align(h(1:6));
irf_zoom(h(1:6),'x',Tint);
set(h(1:6),'fontsize',12);