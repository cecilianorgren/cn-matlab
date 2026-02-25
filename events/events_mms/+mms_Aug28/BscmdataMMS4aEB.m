tsl = 729600;

FID = fopen('mms4_scb_20150828_125300t.dat');
Bscmt = fread(FID,tsl,'int64');
Bscmt = int64(Bscmt);
Bscmt =  EpochTT(Bscmt);

FID = fopen('mms4_scb_20150828_125300Bxyz.dat');
B = fread(FID,tsl*3,'float');
B = [B(1:tsl) B((1+tsl):2*tsl) B((1+2*tsl):3*tsl)];

Bscm = TSeries(Bscmt,B,'to',1);


%load brst mode interval
tmpDataObj = dataobj('data/mms4_edp_brst_ql_dce2d_20150828125314_v0.2.0.cdf');
Exyz = get_variable(tmpDataObj,'mms4_edp_dce_xyz_dsl');
Exyz = mms.variable2ts(Exyz);
Exyz.data(find(abs(Exyz.data) > 200)) = NaN;

tmpDataObj = dataobj('data/mms4_dfg_srvy_ql_20150828_v0.0.3.cdf');
Bxyz = get_variable(tmpDataObj,'mms4_dfg_srvy_dmpa');
Bxyz = mms.variable2ts(Bxyz);

tlimit = irf.tint(Bscm.time.start.utc,Bscm.time.stop.utc);
Exyz = Exyz.tlim(tlimit);
Bscm = Bscm.tlim(tlimit);
tlimitl = tlimit+[-60 60];
Bxyz = Bxyz.tlim(tlimitl);
Bxyz = Bxyz.resample(Exyz);
ecfreq = (1.6e-19)*Bxyz.abs.data*1e-9/(9.1e-31*2*pi);
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = TSeries(Bxyz.time,ecfreq);
ecfreq01 = TSeries(Bxyz.time,ecfreq01);
ecfreq05 = TSeries(Bxyz.time,ecfreq05);
Bscm = Bscm.resample(Exyz);

SCpos = [0 1 0];

Bmag = Bxyz.abs.data;
Rpar = Bxyz.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos);
Rmag   = irf_abs(Rperpy,1);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag   = irf_abs(Rperpx,1);
Rperpx = Rperpx./[Rmag Rmag Rmag];

Epar = dot(Rpar,Exyz.data,2);
Eperp = dot(Rperpx,Exyz.data,2);
Eperp2 = dot(Rperpy,Exyz.data,2);

Bpar = dot(Rpar,Bscm.data,2);
Bperp = dot(Rperpx,Bscm.data,2);
Bperp2 = dot(Rperpy,Bscm.data,2);

Efac = TSeries(Exyz.time,[Eperp Eperp2 Epar],'to',1);
Bfac = TSeries(Bscm.time,[Bperp Bperp2 Bpar],'to',1);

h=irf_plot(3,'newfigure'); 

h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bfac');
irf_plot(h(2),Bfac);
ylabel(h(2),'B_{FAC} (nT)','Interpreter','tex');
irf_legend(h(2),{'B_{\perp 1}','B_{\perp 2}','B_{||}'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Efac');
irf_plot(h(3),Efac);
ylabel(h(3),'E_{FAC} (nT)','Interpreter','tex');
irf_legend(h(3),{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.5 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')


irf_zoom(h(1:3),'x',tlimit);
set(h(1:3),'fontsize',12);

irf_plot_axis_align(h(1:3));
title(h(1),'MMS 4 - Overview 1');

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','BscmoverviewMMS3a.png');
