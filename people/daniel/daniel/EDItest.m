ic = 4;

interval = '20170706135303';
yyyy = interval(1:4);
mm = interval(5:6);
dd = interval(7:8);
starttime = interval(9:end);

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/fgm/brst/l2/' yyyy '/' mm '/' dd '/mms?_fgm_brst_l2_' interval '_v5.92.0.cdf'');'];
c_eval(tmp,ic);
c_eval('Bxyz = mms.variable2ts(get_variable(tmpDataObj,''mms?_fgm_b_gse_brst_l2''));',ic);
c_eval('Bdmpa = mms.variable2ts(get_variable(tmpDataObj,''mms?_fgm_b_dmpa_brst_l2''));',ic);

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/edp/brst/l2/dce/' yyyy '/' mm '/' dd '/mms?_edp_brst_l2_dce_' interval '_v2.2.0.cdf'');'];
c_eval(tmp,ic);
c_eval('Exyz = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_dce_gse_brst_l2''));',ic);

Efac = irf_convert_fac(Exyz, Bxyz, [1, 0, 0]);
[Epar,Eperp]=irf_dec_parperp(Bxyz,Exyz);
Epp = irf.ts_scalar(Epar.time,[Eperp.data Epar.data]);

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/fpi/brst/l2/des-moms/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-moms_' interval '_v3.2.0.cdf'');'];
c_eval(tmp,ic);
c_eval('ne = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_numberdensity_brst''));',ic);
c_eval('Ve = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_bulkv_gse_brst''));',ic);
[Vpar,Vperp]=irf_dec_parperp(Bxyz,Ve);
Vepp = irf.ts_scalar(Vpar.time,[Vperp.data Vpar.data]);

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/edi/brst/l2/amb-pm2/' yyyy '/' mm '/' dd '/mms?_edi_brst_l2_amb-pm2_' interval '_v3.2.0.cdf'');'];
c_eval(tmp,ic);
c_eval('flux0a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_0_brst_l2''));',ic);
c_eval('flux0b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_0_brst_l2''));',ic);

c_eval('flux180a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_180_brst_l2''));',ic);
c_eval('flux180b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_180_brst_l2''));',ic);

flux0180 = irf.ts_scalar(flux180a.time,[(flux0a.data+flux0b.data)/2 (flux180a.data+flux180b.data)/2]);

c_eval('energy1 = get_variable(tmpDataObj,''mms?_edi_energy_gdu1_brst_l2'');',ic);
c_eval('energy2 = get_variable(tmpDataObj,''mms?_edi_energy_gdu2_brst_l2'');',ic);

EDIen = single(energy1.data(1));

EDIdiff = flux0180.data(:,2) - flux0180.data(:,1);
EDIdiff = irf.ts_scalar(flux0180.time,EDIdiff);

EDIdiff180 = EDIdiff.data;
EDIdiff0 = EDIdiff.data;
EDIdiff180(EDIdiff180 < 0) = NaN;
EDIdiff0(EDIdiff0 > 0) = NaN;
EDIdiff = TSeries(EDIdiff.time,[EDIdiff.data EDIdiff180 EDIdiff0]);

%irf_plot({Bxyz,Epp,ne,Vepp,flux0180})

%%

Tint = irf.tint('2017-07-06T13:54:00.00Z/2017-07-06T13:54:10.00Z');

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/edp/brst/l2/scpot/' yyyy '/' mm '/' dd '/mms?_edp_brst_l2_scpot_' interval '_v2.4.0.cdf'');'];
c_eval(tmp,ic);
c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_scpot_brst_l2''));',ic);

tmp = ['tmpDataObj = dataobj(''/Volumes/Nexus/data/mms?/fpi/brst/l2/des-dist/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-dist_' interval '_v3.2.0.cdf'');'];
c_eval(tmp,ic);
ePDist = mms.make_pdist(tmpDataObj);

ePDist = ePDist.tlim(Tint);
Bdmpa = Bdmpa.tlim(Tint);
Bxyz = Bxyz.tlim(Tint);
Vepp = Vepp.tlim(Tint);
ne = ne.tlim(Tint);
Epar = Epar.tlim(Tint);

ePDistomni = ePDist.omni.deflux;
ePDistpitch = ePDist.pitchangles(Bdmpa,4).deflux;

PSDparapar = ePDistpitch.data(:,:,1)./ePDistpitch.data(:,:,end);
specparapar =struct('t',ePDistpitch.time.epochUnix);
specparapar.p = PSDparapar;
specparapar.p_label={'log_{10}(f_{||+}/f_{||-})'};
specparapar.f_label={''};
specparapar.f = ePDistomni.depend{1,1};

%%

ediline = irf.ts_scalar(Tint,[EDIen EDIen]');

h=irf_plot(7,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.86;
ywidth = 0.13;
set(h(1),'position',[0.10 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.10 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.10 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.10 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.10 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.10 0.97-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.10 0.97-7*ywidth xwidth ywidth]);

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),{'B (nT)'},'Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.98 0.12])
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',12)

h(2)=irf_panel('edist');
irf_spectrogram(h(2),ePDistomni.specrec,'log');
hold(h(2),'on');
irf_plot(h(2),SCpot,'k')
irf_plot(h(2),ediline,'m');
hold(h(2),'off');
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(2),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(3)=irf_panel('edistparapar');
irf_spectrogram(h(3),specparapar,'log','donotfitcolorbarlabel');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
hold(h(3),'on');
irf_plot(h(3),SCpot,'k');
irf_plot(h(3),ediline,'m');
hold(h(3),'off');
set(h(3),'yscale','log');
caxis(h(3),[-2 2])
set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(3),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(4)=irf_panel('EDI');
irf_plot(h(4),EDIdiff);
ylabel(h(4),{'flux 180^{o}-0^{o}','(cm^2 s^{-1})'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.90],'color','k','fontsize',12)

h(5)=irf_panel('Vepp');
irf_plot(h(5),Vepp);
ylabel(h(5),{'V (km s^{-1})'},'Interpreter','tex');
irf_legend(h(5),{'V_{x,\perp}','V_{y,\perp}','V_{z,\perp}','V_{||}'},[0.98 0.12])
irf_legend(h(5),'(e)',[0.99 0.94],'color','k','fontsize',12)

h(6)=irf_panel('ne');
irf_plot(h(6),ne);
ylabel(h(6),{'n_e (cm^{-3})'},'Interpreter','tex');
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',12)

h(7)=irf_panel('Epar');
irf_plot(h(7),Epar);
ylabel(h(7),{'E_{||} (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(7),'(g)',[0.99 0.94],'color','k','fontsize',12)

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

colormap(h(2),'jet');
colormap(h(3),bgrcmap);

c_eval('title(h(1),''MMS?'');',ic)

irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',Tint);
set(h(1:7),'fontsize',12);
