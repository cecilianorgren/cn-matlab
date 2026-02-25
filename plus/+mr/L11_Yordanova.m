tint_caa = irf.tint('2003-04-16T20:30:00.00Z/2003-04-16T21:30:00.00Z');
time = irf_time('2003-04-16T21:18:16.00Z','utc>epochtt');

caa_download(tint_caa,'C?_CP_FGM_FULL');
%%
c_eval('gseB?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''ts'');');


c_eval('gseB?.userData.LABLAXIS = ''B''; gseB?.name=''B (nT)'';',1:4);
%h = irf_plot({gseB1,gseB2,gseB3,gseB4});
%c_eval('irf_pl_mark(h(?),time,''g'')',1:4)

c_eval('dB?_1 = irf.ts_vec_xyz(gseB?.time(2:end),diff(gseB?.data,1)); dB?_1.name = ''dB (nT)'';',1:4)
c_eval('dB?_5 = irf.ts_vec_xyz(gseB?.time(3:1:end-2),gseB?.data(5:end,:)-gseB?.data(1:end-4,:)); dB?_5.name = ''dB (nT)'';',1:4)
c_eval('dB?_11 = irf.ts_vec_xyz(gseB?.time(6:1:end-5),gseB?.data(11:end,:)-gseB?.data(1:end-10,:)); dB?_11.name = ''dB (nT)'';',1:4)
c_eval('dB?_21 = irf.ts_vec_xyz(gseB?.time(11:1:end-10),gseB?.data(21:end,:)-gseB?.data(1:end-20,:)); dB?_21.name = ''dB (nT)'';',1:4)

%c_eval('thB?_1 = irf.ts_scalar(gseB?.time(2:end),(sum(gseB?.data(2:end,:).*gseB?.data(1:end-1,:),2)./(abs(gseB?.data(2:end,:)).*abs(gseB?.data(1:end-1,:)))); thB?_1.name = ''thB (nT)'';',1:4)
%c_eval('thB?_5 = irf.ts_scalar(gseB?.time(3:1:end-2),(sum(gseB?.data(5:end,:).*gseB?.data(1:end-4,:),2)./(abs(gseB?.data(5:end,:)).*abs(gseB?.data(1:end-4,:)))); thB?_5.name = ''thB (nT)'';',1:4)
%c_eval('thB?_11 = irf.ts_scalar(gseB?.time(6:1:end-5),(sum(gseB?.data(11:end,:).*gseB?.data(1:end-10,:),2)./(abs(gseB?.data(11:end,:)).*abs(gseB?.data(1:end-10,:)))); thB?_11.name = ''thB (nT)'';',1:4)

thB1_11 = irf.ts_scalar(gseB1.time(6:1:end-5),acosd(sum(gseB1.data(11:end,:).*gseB1.data(1:end-10,:),2)./((gseB1.abs.data(11:end,:)).*(gseB1.abs.data(1:end-10,:))))); 
c_eval('thB?_11 = irf.ts_scalar(gseB?.time(6:1:end-5),acosd(sum(gseB?.data(11:end,:).*gseB?.data(1:end-10,:),2)./((gseB?.abs.data(11:end,:)).*(gseB?.abs.data(1:end-10,:)))));',1:4)
c_eval('thB?_21 = irf.ts_scalar(gseB?.time(11:1:end-10),acosd(sum(gseB?.data(21:end,:).*gseB?.data(1:end-20,:),2)./((gseB?.abs.data(21:end,:)).*(gseB?.abs.data(1:end-20,:)))));',1:4)

c_eval('pviB?_1 = dB?_1.abs/std(dB?_1.abs.data); pviB?_1.name = ''PVI'';',1:4)
c_eval('pviB?_5 = dB?_5.abs/std(dB?_5.abs.data); pviB?_5.name = ''PVI'';',1:4)
c_eval('pviB?_11 = dB?_11.abs/std(dB?_11.abs.data); pviB?_11.name = ''PVI'';',1:4)
c_eval('pviB?_21 = dB?_11.abs/std(dB?_21.abs.data); pviB?_21.name = ''PVI'';',1:4)
%%
h = irf_plot(3);
hca = irf_panel('B1');
irf_plot(hca,gseB1)
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.02, 0.8]);

hca = irf_panel('PVI_11 B1');
irf_plot(hca,pviB1_11)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B1');
irf_plot(hca,thB1_11)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h,'x',tint_caa)
irf_plot_axis_align

hca = irf_panel('th_11 B1'); hca.YLim = [0 180];

for ip = 1:3
    h(ip).FontSize = 14;
end
%% C3, 11
h = irf_plot(3);
hca = irf_panel('B3');
irf_plot(hca,gseB3)
hca.YLabel.String = 'B (nT)';

hca = irf_panel('PVI_11 B3');
irf_plot(hca,pviB3_11)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B3');
irf_plot(hca,thB3_11)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h,'x',tint_caa)
irf_plot_axis_align

hca = irf_panel('th_11 B3'); hca.YLim = [0 180];

for ip = 1:3
    h(ip).FontSize = 14;
end
%% C3, 21
h = irf_plot(3);
hca = irf_panel('B3');
irf_plot(hca,gseB3)
hca.YLabel.String = 'B (nT)';

hca = irf_panel('PVI_11 B3');
irf_plot(hca,pviB3_21)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B3');
irf_plot(hca,thB3_21)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h,'x',tint_caa)
irf_plot_axis_align

hca = irf_panel('th_11 B3'); hca.YLim = [0 180];

for ip = 1:3
    h(ip).FontSize = 14;
end
%%

scatter(thB1_11.data,pviB1_11.data,'.')
%%
hca = irf_panel('B2');
irf_plot(hca,gseB2)
hca = irf_panel('PVI B2');
irf_plot(hca,pviB2)
hca = irf_panel('B3');
irf_plot(hca,gseB3)
hca = irf_panel('PVI B3');
irf_plot(hca,pviB3)
hca = irf_panel('B4');
irf_plot(hca,gseB4)
hca = irf_panel('PVI B4');
irf_plot(hca,pviB4)

%% Zoom in on large pvi
tint_pvi = irf.tint('2003-04-16T21:09:00.00Z/2003-04-16T21:09:10.00Z');
tint_mva = irf.tint('2003-04-16T21:09:02.50Z/2003-04-16T21:09:06.90Z');

[mvaB1,l,v]=irf_minvar(gseB1.tlim(tint_mva));

h = irf_plot(3);
hca = irf_panel('B1');
irf_plot(hca,gseB1)
hca.YLabel.String = 'B (nT)';

hca = irf_panel('mvaB1');
irf_plot(hca,{mvaB1.x,mvaB1.y,mvaB1.z,mvaB1.abs},'comp')
hca.YLabel.String = 'B (nT)';

hca = irf_panel('PVI_11 B1');
irf_plot(hca,pviB1_11)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B1');
irf_plot(hca,thB1_11)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h,'x',tint_mva)
irf_plot_axis_align

hca = irf_panel('th_11 B1'); hca.YLim = [0 180];

for ip = 1:3
    h(ip).FontSize = 14;
end

%% Zoom in on large pvi, mva
tint_pvi = irf.tint('2003-04-16T21:09:00.00Z/2003-04-16T21:09:10.00Z');
tint_mva = irf.tint('2003-04-16T21:09:02.50Z/2003-04-16T21:09:06.90Z');

[mvaB1,l,v]=irf_minvar(gseB1.tlim(tint_mva));

[h1,h2] = initialize_combined_plot(3,2,1,0.6,'vertical');
%%
hca = irf_panel('mvaB1');
irf_plot(hca,{mvaB1.x,mvaB1.y,mvaB1.z,mvaB1.abs},'comp')
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.02, 0.8]);

hca = irf_panel('PVI_11 B1');
irf_plot(hca,pviB1_11)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B1');
irf_plot(hca,thB1_11)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h1,'x',tint_mva)
irf_plot_axis_align

hca = irf_panel('th_11 B1'); hca.YLim = [0 180];

for ip = 1:3
   h1(ip).FontSize = 14;
end

plot(h2(1),mvaB1.x.data,mvaB1.z.data);
h2(1).XLabel.String = 'max';
h2(1).YLabel.String = 'min';
h2(1).Title.String = sprintf('l_1/l_2 = %.2f, l_2/l_3 = %.2f',l(1)/l(2),l(2)/l(3));

plot(h2(2),mvaB1.x.data,mvaB1.y.data);
h2(2).XLabel.String = 'max';
h2(2).YLabel.String = 'inter';


for ip = 1:2
    axis(h2(ip),'equal')
    h2(ip).XGrid = 'on';
    h2(ip).YGrid = 'on';
end

%%
hca = irf_panel('mvaB1');
irf_plot(hca,{mvaB1.x,mvaB1.y,mvaB1.z,mvaB1.abs},'comp')
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.02, 0.8]);

hca = irf_panel('PVI_11 B1');
irf_plot(hca,pviB1_11)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_11 B1');
irf_plot(hca,thB1_11)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h1,'x',tint_mva)
irf_plot_axis_align

hca = irf_panel('th_11 B1'); hca.YLim = [0 180];

for ip = 1:3
   h1(ip).FontSize = 14;
end

plot(h2(1),mvaB1.x.data,mvaB1.z.data);
h2(1).XLabel.String = 'max';
h2(1).YLabel.String = 'min';
h2(1).Title.String = sprintf('l_1/l_2 = %.2f, l_2/l_3 = %.2f',l(1)/l(2),l(2)/l(3));

plot(h2(2),mvaB1.x.data,mvaB1.y.data);
h2(2).XLabel.String = 'max';
h2(2).YLabel.String = 'inter';


for ip = 1:2
    axis(h2(ip),'equal')
    h2(ip).XGrid = 'on';
    h2(ip).YGrid = 'on';
end
    
%% Zoom in on large pvi, mva
ic = 2;
tint_pvi = irf.tint('2003-04-16T21:09:00.00Z/2003-04-16T21:09:10.00Z');
tint_mva = irf.tint('2003-04-16T21:18:14.50Z/2003-04-16T21:18:17.50Z');

c_eval('[mvaB?,l,v]=irf_minvar(gseB?.tlim(tint_mva));',ic)

[h1,h2] = initialize_combined_plot(3,2,1,0.6,'vertical');


hca = irf_panel('mvaB');
c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z,mvaB?.abs},''comp'')',ic)
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.02, 0.8]);

hca = irf_panel('PVI_21');
c_eval('irf_plot(hca,pviB?_21)',ic)
hca.YLabel.String = 'PVI';

hca = irf_panel('th_21');
c_eval('irf_plot(hca,thB?_21)',ic)
hca.YLabel.String = '\theta_B (\circ)';
ylabel(hca,'\theta_B (\circ)','interpreter','tex')

irf_zoom(h1,'x',tint_mva)
irf_plot_axis_align

hca = irf_panel('th_11'); hca.YLim = [0 180];

for ip = 1:3
   h1(ip).FontSize = 14;
end

c_eval('plot(h2(1),mvaB?.x.data,mvaB?.z.data);',ic)
h2(1).XLabel.String = 'max';
h2(1).YLabel.String = 'min';
h2(1).Title.String = sprintf('l_1/l_2 = %.2f, l_2/l_3 = %.2f',l(1)/l(2),l(2)/l(3));

c_eval('plot(h2(2),mvaB?.x.data,mvaB?.y.data);',ic)
h2(2).XLabel.String = 'max';
h2(2).YLabel.String = 'inter';


for ip = 1:2
    axis(h2(ip),'equal')
    h2(ip).XGrid = 'on';
    h2(ip).YGrid = 'on';
end   
    