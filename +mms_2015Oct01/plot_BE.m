%tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T13:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
clustercolors = [0 0 0;1 0 0; 0 1 0;0 0 1];

h = irf_plot(7);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 1;


hca = irf_panel('B fg');
c_eval('irf_plot(hca,{dmpaB?.x,dmpaB?.y,dmpaB?.z,dmpaB?.abs},''comp'');',ic)
hca.YLabel.String = {'B_{fg,DMPA}','[nT]'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('B sc');
c_eval('irf_plot(hca,{B?sc.x,B?sc.y,B?sc.z,B?sc.abs},''comp'');',ic)
hca.YLabel.String = {'B_{sc}','[nT]'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('E');
c_eval('irf_plot(hca,{dslE?.x,dslE?.y,dslE?.z,dslE?.abs},''comp'');',ic)
hca.YLabel.String = {'E_{DSL}','[mV/m]'};
irf_legend(hca,{'E_x','E_y','E_z','|E|'},[0.95 0.95]);

hca = irf_panel('E brst');
c_eval('irf_plot(hca,{dslE?brst.x,dslE?brst.y,dslE?brst.z,dslE?brst.abs},''comp'');',ic)
hca.YLabel.String = {'E_{DSL}','[mV/m]'};
irf_legend(hca,{'E_x','E_y','E_z','|E|'},[0.95 0.95]);


hca = irf_panel('-scPot');
c_eval('irf_plot(hca,P?*(-1));',ic)
hca.YLabel.String = {'-scPot','[V]'};
hca.YScale = 'lin';

hca = irf_panel('-scPot brst');
c_eval('irf_plot(hca,P?brst*(-1));',ic)
hca.YLabel.String = {'-scPot','[V]'};
hca.YScale = 'lin';

if 1
    hca = irf_panel('J');
    irf_plot(hca,j);
    hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);
end
if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end

tintOv = irf.tint('2015-10-01T06:53:00.00Z/2015-10-01T06:56:10.00Z');
tint1 = irf.tint('2015-10-01T06:55:33.000Z/2015-10-01T06:55:35.000Z'); % sharp crossing where fields are similar
tint2 = irf.tint('2015-10-01T06:53:32.000Z/2015-10-01T06:53:48.000Z'); % the one with so different fields

irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')
irf_zoom(h,'x',tintOv)
irf_zoom(h,'y')
irf_plot_axis_align
%%

h = irf_plot(7);

hca = irf_panel('Bx');
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y,DMPA}','[nT]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{dslE1brst.x.tlim(tint),dslE2brst.x.tlim(tint),dslE3brst.x.tlim(tint),dslE4brst.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey');
irf_plot(hca,{dslE1brst.y.tlim(tint),dslE2brst.y.tlim(tint),dslE3brst.y.tlim(tint),dslE4brst.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y,DSL}','[mV/m]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('-Vsc');
%c_eval('mP? = P?; mP?.data = -mP?.data;',ic)
irf_plot(hca,{P1brst*-1,P2brst*-1,P3brst*-1,P4brst*-1},'comp');
hca.YLabel.String = {'-V_{SC}','[V]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
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

