tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T13:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
clustercolors = [0 0 0;1 0 0; 0 1 0;0 0 1];

h = irf_plot(8);
for ii = 1:8
    set(h(ii),'ColorOrder',clustercolors)
end


hca = irf_panel('Bx');
irf_plot(hca,{dfg1.x.tlim(tint),dfg2.x.tlim(tint),dfg3.x.tlim(tint),dfg4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('By');
irf_plot(hca,{dfg1.y.tlim(tint),dfg2.y.tlim(tint),dfg3.y.tlim(tint),dfg4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('Bz');
irf_plot(hca,{dfg1.z.tlim(tint),dfg2.z.tlim(tint),dfg3.z.tlim(tint),dfg4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('Ex');
irf_plot(hca,{edp1.x.tlim(tint),edp2.x.tlim(tint),edp3.x.tlim(tint),edp4.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('Ey');
irf_plot(hca,{edp1.y.tlim(tint),edp2.y.tlim(tint),edp3.y.tlim(tint),edp4.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('-Vsc');
%c_eval('mP? = P?; mP?.data = -mP?.data;',sc)
irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
hca.YLabel.String = {'-V_{SC}','[V]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95],'color','cluster');

hca = irf_panel('J');
irf_plot(hca,j*1e9);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{DMPA}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

hca = irf_panel('Jfac');
irf_plot(hca,jfac*1e9);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);



irf_zoom(h,'x',tint)
irf_zoom(h,'y')

%%
tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T10:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
clustercolors = [0 0 0;0 0 1; 1 0 0;0.3 0.3 0.3];

h = irf_plot(8);



hca = irf_panel('Bx');
irf_plot(hca,{dfg1.x.tlim(tint),dfg2.x.tlim(tint),dfg3.x.tlim(tint),dfg4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{dfg1.y.tlim(tint),dfg2.y.tlim(tint),dfg3.y.tlim(tint),dfg4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{dfg1.z.tlim(tint),dfg2.z.tlim(tint),dfg3.z.tlim(tint),dfg4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{edp1.x.tlim(tint),edp2.x.tlim(tint),edp3.x.tlim(tint),edp4.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey');

irf_plot(hca,{edp1.y.tlim(tint),edp2.y.tlim(tint),edp3.y.tlim(tint),edp4.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('-Vsc');
c_eval('mP? = P?; mP?.data = -mP?.data;',sc)
irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
hca.YLabel.String = {'-V_{SC}','[V]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j*1e9);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{DMPA}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

hca = irf_panel('Jfac');
irf_plot(hca,jfac*1e9);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);



irf_zoom(h,'x',tint)
irf_zoom(h,'y')

