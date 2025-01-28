tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');

c_eval('tmpDataObj? = dataobj(''data/mms?_dfg_brst_l2pre_20151016103254_v2.5.0.cdf'');',[1:4]);
c_eval('Bxyz? = get_variable(tmpDataObj?,''mms?_dfg_brst_l2pre_gse'');',[1:4]);
c_eval('Bxyz? = mms.variable2ts(Bxyz?);',[1:4]);
c_eval('Bxyz? = TSeries(Bxyz?.time,Bxyz?.data(:,[1:3]),''to'',1);',[1:4]);
c_eval('Bxyz? = Bxyz?.tlim(tint);',[1:4]);
c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',[2:4]);
Bxyzav = (Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4;
Bxyzav = TSeries(Bxyz1.time,Bxyzav,'to',1);
%tint = irf.tint(Bxyz1.time.start.utc,Bxyz1.time.stop.utc);
if 1,
c_eval('Rxyz? = get_variable(tmpDataObj?,''mms?_pos_gse'');',[1:4]);
c_eval('Rxyz? = mms.variable2ts(Rxyz?);',[1:4]);
c_eval('Rxyz? = TSeries(Rxyz?.time,Rxyz?.data(:,[1:3]),''to'',1);',[1:4]);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',[1:4]);
end
if 0,
R  = mms.get_data('R_gse',tint);
c_eval('Rxyz? = TSeries(R.time,R.gseR?,''to'',1);',[1:4]);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',[1:4]);
end

c_eval('tmpDataObj? = dataobj(''data/mms?_edp_fast_ql_dce2d_20151016000000_v0.4.3.cdf'');',[1:4]);
c_eval('Exyz? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_edp_dce_xyz_dsl''));',[1:4]);
c_eval('Exyz? = Exyz?.tlim(tint);',[1:4]);
c_eval('Exyz? = Exyz?.resample(Exyz1);',[2:4]);
Exyzav = (Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)/4;
Exyzav = TSeries(Exyz1.time,Exyzav,'to',1);
Exyzav = Exyzav+[-2 0 0]; % Ex offset?

c_eval('tmpDataObj? = dataobj(''data/mms?_fpi_brst_l1b_des-moms_20151016103300_v1.0.0.cdf'');',[1:4]);
c_eval('ne? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_numberDensity''));',[1:4]);
c_eval('dtmoments? = ne?.time(2)-ne?.time(1);',[1:4]);
c_eval('fsmoments? = 1/dtmoments?;',[1:4]);
c_eval('ne? = irf_filt(ne?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('ne? = ne?.tlim(tint);',[1:4]);
c_eval('ne? = ne?.resample(Bxyz1);',[1:4]);
neav = TSeries(ne1.time,(ne1.data+ne2.data+ne3.data+ne4.data)/4);
c_eval('TeX? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_TempXX''));',[1:4]);
c_eval('TeY? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_TempYY''));',[1:4]);
c_eval('TeZ? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_TempZZ''));',[1:4]);

c_eval('PeXX? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresXX''));',[1:4]);
c_eval('PeXY? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresXY''));',[1:4]);
c_eval('PeXZ? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresXZ''));',[1:4]);
c_eval('PeYY? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresYY''));',[1:4]);
c_eval('PeYZ? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresYZ''));',[1:4]);
c_eval('PeZZ? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_PresZZ''));',[1:4]);
c_eval('PeXX? = irf_filt(PeXX?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeXY? = irf_filt(PeXY?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeXZ? = irf_filt(PeXZ?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeYY? = irf_filt(PeYY?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeYZ? = irf_filt(PeYZ?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeZZ? = irf_filt(PeZZ?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('PeXX?.data = PeXX?.data*1e-9;',[1:4]); % Conversion to SI units for v1.0.0 data
c_eval('PeXY?.data = PeXY?.data*1e-9;',[1:4]);
c_eval('PeXZ?.data = PeXZ?.data*1e-9;',[1:4]);
c_eval('PeYY?.data = PeYY?.data*1e-9;',[1:4]);
c_eval('PeYZ?.data = PeYZ?.data*1e-9;',[1:4]);
c_eval('PeZZ?.data = PeZZ?.data*1e-9;',[1:4]);
c_eval('PeXX? = PeXX?.resample(Bxyz1);',[1:4]);
c_eval('PeXY? = PeXY?.resample(Bxyz1);',[1:4]);
c_eval('PeXZ? = PeXZ?.resample(Bxyz1);',[1:4]);
c_eval('PeYY? = PeYY?.resample(Bxyz1);',[1:4]);
c_eval('PeYZ? = PeYZ?.resample(Bxyz1);',[1:4]);
c_eval('PeZZ? = PeZZ?.resample(Bxyz1);',[1:4]);

c_eval('VeX? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_bulkX''));',[1:4]);
c_eval('VeY? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_bulkY''));',[1:4]);
c_eval('VeZ? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_bulkZ''));',[1:4]);
c_eval('Ve? = TSeries(VeX?.time,[VeX?.data VeY?.data VeZ?.data],''to'',1);',[1:4]);
%Veav = TSeries(VeX1.time,(Ve1.data+Ve2.data+Ve3.data+Ve4.data)/4,'to',1);
%Veav = Veav.resample(Bxyz1);
%EvxBe = cross(Veav,Bxyzav);
%EvxBe.data = -EvxBe.data*1e-3;
c_eval('Ve? = Ve?.resample(Bxyz1);',[1:4]);
c_eval('EvxBe? = cross(Ve?,Bxyz?);',[1:4]);
c_eval('EvxBe?.data = -EvxBe?.data*1e-3;',[1:4]);
EvxBe = TSeries(EvxBe1.time,(EvxBe1.data+EvxBe2.data+EvxBe3.data+EvxBe4.data)/4,'to',1);


EPeXX = c_4_grad('Rxyz?','PeXX?','grad');
EPeXY = c_4_grad('Rxyz?','PeXY?','grad');
EPeXZ = c_4_grad('Rxyz?','PeXZ?','grad');
EPeYY = c_4_grad('Rxyz?','PeYY?','grad');
EPeYZ = c_4_grad('Rxyz?','PeYZ?','grad');
EPeZZ = c_4_grad('Rxyz?','PeZZ?','grad');
EPeX = -(EPeXX.data(:,1)+EPeXY.data(:,2)+EPeXZ.data(:,3))./(neav.data*1e6*1.6e-19);
EPeY = -(EPeXY.data(:,1)+EPeYY.data(:,2)+EPeYZ.data(:,3))./(neav.data*1e6*1.6e-19);
EPeZ = -(EPeXZ.data(:,1)+EPeYZ.data(:,2)+EPeZZ.data(:,3))./(neav.data*1e6*1.6e-19);
EPe = TSeries(EPeXX.time,[EPeX EPeY EPeZ],'to',1);
EPefac = irf_convert_fac(EPe,Bxyzav,Rxyz1);

[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');
j.data(:,2) = j.data(:,2)-70e-9; % Remove some offsets?
j.data(:,1) = j.data(:,1)-25e-9;
jxB = cross(j,Bxyzav);
jxB.data = jxB.data*1e-9;
jxB.data = jxB.data./[neav.data neav.data neav.data]; 
jxB.data = jxB.data/1.6e-19/1000; %Convert to (mV/m)

c_eval('tmpDataObj? = dataobj(''data/mms?_fpi_brst_l1b_dis-moms_20151016103300_v1.0.0.cdf'');',[1:4]);
c_eval('ni? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_numberDensity''));',[1:4]);
c_eval('dtmoments? = ni?.time(2)-ni?.time(1);',[1:4]);
c_eval('fsmoments? = 1/dtmoments?;',[1:4]);
c_eval('ni? = irf_filt(ni?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('Vx? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_bulkX''));',[1:4]);
c_eval('Vy? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_bulkY''));',[1:4]);
c_eval('Vz? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_bulkZ''));',[1:4]);
c_eval('Vi? = TSeries(Vx?.time,[Vx?.data Vy?.data Vz?.data],''to'',1);',[1:4]);
c_eval('Vi? = irf_filt(Vi?,0,fsmoments?/4,fsmoments?,3);',[1:4]);
c_eval('Vi? = Vi?.resample(Bxyzav);',[1:4]);
%Viav = TSeries(Vi1.time,(Vi1.data+Vi2.data+Vi3.data+Vi4.data)/4,'to',1);
%EvxBi = cross(Viav,Bxyzav);
%EvxBi.data = -EvxBi.data*1e-3;
c_eval('EvxBi? = cross(Vi?,Bxyz?);',[1:4]);
c_eval('EvxBi?.data = -EvxBi?.data*1e-3;',[1:4]);
EvxBi = TSeries(EvxBi1.time,(EvxBi1.data+EvxBi2.data+EvxBi3.data+EvxBi4.data)/4,'to',1);

Etot = TSeries(EvxBi.time,EvxBi.data+jxB.data+EPe.data,'to',1);

h = irf_plot(7,'newfigure');

hca = irf_panel('BMMS1');
irf_plot(hca,Bxyz1);
ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca = irf_panel('EMMS1');
irf_plot(hca,Exyzav);
ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('jxBfield');
irf_plot(hca,jxB);
ylabel(hca,{'J \times B/q_{e} n_{e}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca = irf_panel('EvxB');
irf_plot(hca,EvxBi);
ylabel(hca,{'-V_{i} \times B','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(d)',[0.99 0.98],'color','k')

hca = irf_panel('Ep');
irf_plot(hca,EPe);
ylabel(hca,{'- \nabla \cdot P_{e}/q_{e} n_{e}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(e)',[0.99 0.98],'color','k')

hca = irf_panel('Etot');
irf_plot(hca,Etot);
ylabel(hca,{'Etot','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(f)',[0.99 0.98],'color','k')

hca = irf_panel('EvxBe');
irf_plot(hca,EvxBe);
ylabel(hca,{'-V_{e} \times B','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(g)',[0.99 0.98],'color','k')

title(h(1),'MMS - 4 Spacecraft average')

irf_plot_axis_align(1,h(1:7))
irf_zoom(h(1:7),'x',tint);