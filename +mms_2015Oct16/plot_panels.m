% Example panels that can be copied and pasted
%% Initialize figure

h = irf_plot(7); 
isub = 1;
hca = irf_panel('Example'); isub = isub + 1;

%% Finishing touches to figure
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

%% dmpaB?brst 1sc (ic)
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);
%% mvaB? 4sc
hca = irf_panel('B_L');
irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_M');
irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{M}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_N');
irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{N}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
%% dmpaB? 4sc
hca = irf_panel('Bx');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
%% dslE?brst 1sc (ic)
hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[0.95 0.95]);

%%
tint = irf.tint('2015-10-16T10:33:22.000Z/2015-10-16T10:33:38.000Z');
h = irf_plot(6);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 2;
%c_eval('mvaEdb?.data(abs(ang?.data)<20,3)=NaN;')


hca = irf_panel('B_L');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('B_M');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{M}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('B_N');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{N}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('Babs');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
hca.YLabel.String = {'|B|','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('J');
set(hca,'ColorOrder',mms_colors('xyza'))
irf_plot(hca,{mvaJ.x,mvaJ.y,mvaJ.z},'comp');
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'J_L','J_M','J_N'},[0.95 0.95]);


hca = irf_panel('mvaE?');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic);
hca.YLabel.String = {'E','[mV/m]'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_L','E_M','E_N'},[0.95 0.95]);
irf_legend(hca,{irf_ssub('mms?',ic)},[0.05 0.95],'color',mms_colors(irf_ssub('?',ic)));



if 0 % spacecraft potential 
    hca = irf_panel('-scPot');
    irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
    hca.YLabel.String = {'-scPot','[V]'};
    irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
    hca.YScale = 'lin';
end



if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end


