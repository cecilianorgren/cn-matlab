%tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T13:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
clustercolors = [0 0 0;1 0 0; 0 1 0;0 0 1];

h = irf_plot(8);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 1;


hca = irf_panel('Bx');
irf_plot(hca,{Bxyz1.x.tlim(tint),Bxyz2.x.tlim(tint),Bxyz3.x.tlim(tint),Bxyz4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{Bxyz1.y.tlim(tint),Bxyz2.y.tlim(tint),Bxyz3.y.tlim(tint),Bxyz4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{Bxyz1.z.tlim(tint),Bxyz2.z.tlim(tint),Bxyz3.z.tlim(tint),Bxyz4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Babs');
irf_plot(hca,{Bxyz1.abs.tlim(tint),Bxyz2.abs.tlim(tint),Bxyz3.abs.tlim(tint),Bxyz4.abs.tlim(tint)},'comp');
hca.YLabel.String = {'|B|','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{Exyz1.x.tlim(tint),Exyz2.x.tlim(tint),Exyz3.x.tlim(tint),Exyz4.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey');
irf_plot(hca,{Exyz1.y.tlim(tint),Exyz2.y.tlim(tint),Exyz3.y.tlim(tint),Exyz4.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('-scPot');
irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
hca.YLabel.String = {'-scPot','[V]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
hca.YScale = 'lin';

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
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
hca.YLabel.String = {'E_{y}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('-Vsc');
c_eval('mP? = P?; mP?.data = -mP?.data;',ic)
irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
hca.YLabel.String = {'-V_{SC}','[V]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
%hca.Yscale = 'log';

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')

