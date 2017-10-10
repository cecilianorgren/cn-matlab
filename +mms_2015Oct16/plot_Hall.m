%tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T13:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];

if 0
    %%
    tint = irf.tint('2015-10-16T10:33:00.000Z/2015-10-16T10:36:00.000Z');
    irf_minvar_gui(dmpaB3.tlim(tint))
end
if 1
    % Automatic results made from gui above
    ic = 1:4;
    tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
    c_eval('[out?,l?,v?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
    v = v3;
    % rotate around x-axis 
    turn_angle = 0; % degrees
    tmpv = v;
    v(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
    v(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
     % v is 3x3 matrix:
     %v1: [0.0906 0.0896 0.9918] - first row
     %v2: [0.3510 -0.9349 0.0524] - second row
     %v3: [0.9320 0.3434 -0.1161] - third row   
end
tint = tint_mva;
%mvaB1 = 
c_eval('mvaB?=irf.ts_vec_xyz(dmpaB?.time,[dmpaB?.dot(v(1,:)).data_ dmpaB?.dot(v(2,:)).data_ dmpaB?.dot(v(3,:)).data_]);',ic)
c_eval('mvaE?=irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(v(1,:)).data_ dslE?brst.dot(v(2,:)).data_ dslE?brst.dot(v(3,:)).data_]);')
mvaJ=irf.ts_vec_xyz(j.time,[j.data*v(1,:)' j.data*v(2,:)' j.data*v(3,:)']);
c_eval('[EdB?, ang?] = irf_edb(irf.ts_vec_xy(dslE?brst.time,dslE?brst.data(:,1:2)),dmpaB?); EdB?.units = ''mV/m'';')
c_eval('mvaEdb?=irf.ts_vec_xyz(EdB?.time,[EdB?.data*v(1,:)'' EdB?.data*v(2,:)'' EdB?.data*v(3,:)'']);',ic)
c_eval('mvaEJxB?=irf.ts_vec_xyz(EJxB?.time,[EJxB?.data*v(1,:)'' EJxB?.data*v(2,:)'' EJxB?.data*v(3,:)'']);',ic)

%irf_plot({mvaB1.tlim(tint),mvaB2.tlim(tint),mvaB3.tlim(tint),mvaB4.tlim(tint)},'comp')
%%
h = irf_plot(8);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 2;


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

hca = irf_panel('EJxB');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{mvaEJxB?.x,mvaEJxB?.y,mvaEJxB?.z},''comp'');',ic)
ylabel(hca,{'j\times B/(ne)','[mV/m]'},'interpreter','tex')
%hca.YLabel.String = {'j\times B/(ne)','[mV/m]'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_L','E_M','E_N'},[0.95 0.95]);
irf_legend(hca,{irf_ssub('mms?',ic)},[0.05 0.95]);

hca = irf_panel('mvaE?');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{mvaEdb?.x,mvaEdb?.y,mvaEdb?.z},''comp'');',ic);
hca.YLabel.String = {'E','[mV/m]'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_L','E_M','E_N'},[0.95 0.95]);
irf_legend(hca,{irf_ssub('mms?',ic)},[0.05 0.95]);

hca = irf_panel('x terms');
set(hca,'ColorOrder',mms_colors('xyz'))
c_eval('irf_plot(hca,{dslE?.x,vexB?.x,gradPene?.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'X-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);


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


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

%%
h = irf_plot(7);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 1;


hca = irf_panel('B_L');
irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_M');
irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_N');
irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{N}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('JxB_L');
irf_plot(hca,{mvaEJxB1.x.tlim(tint),mvaEJxB2.x.tlim(tint),mvaEJxB3.x.tlim(tint)},'comp');
hca.YLabel.String = {'JxB_{L}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('JxB_M');
irf_plot(hca,{mvaEJxB1.y.tlim(tint),mvaEJxB2.y.tlim(tint),mvaEJxB3.y.tlim(tint),mvaEJxB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'JxB_{M}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('JxB_N');
irf_plot(hca,{mvaEJxB1.z.tlim(tint),mvaEJxB2.z.tlim(tint),mvaEJxB3.z.tlim(tint),mvaEJxB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'JxB_{N}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('EJxB');
irf_plot(hca,mvaEJxB3);
ylabel(hca,{'j\times B/(ne)','[mV/m]'},'interpreter','tex')
%hca.YLabel.String = {'j\times B/(ne)','[mV/m]'};
irf_legend(hca,{'E_L','E_M','E_N'},[0.95 0.95]);

if 0 % spacecraft potential 
    hca = irf_panel('-scPot');
    irf_plot(hca,{mP1,mP2,mP3,mP4},'comp');
    hca.YLabel.String = {'-scPot','[V]'};
    irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
    hca.YScale = 'lin';
end

hca = irf_panel('J');
irf_plot(hca,mvaJ);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_L','J_M','J_N'},[0.95 0.95]);

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
h = irf_plot(5);
%for ii = 1:8
%    set(h(ii),'ColorOrder',clustercolors)
%end
ic = 1;


hca = irf_panel('B?');
c_eval('irf_plot(hca,mvaB?);',ic);
hca.YLabel.String = {irf_ssub('B_?',ic),'[nT]'};
irf_legend(hca,{'B_L','B_M','B_N'},[0.95 0.95]);


hca = irf_panel('JxB L');
c_eval('irf_plot(hca,{mvaEJxB?.x,mvaEdb?.x},''comp'');',ic);
hca.YLabel.String = {'JxB_{L}','[mV/m]'};
irf_legend(hca,{'JxB','E'},[0.95 0.95]);

hca = irf_panel('JxB M');
c_eval('irf_plot(hca,{mvaEJxB?.y,mvaEdb?.y},''comp'');',ic);
hca.YLabel.String = {'JxB_{M}','[mV/m]'};
irf_legend(hca,{'JxB','E'},[0.95 0.95]);

hca = irf_panel('JxB N');
c_eval('irf_plot(hca,{mvaEJxB?.z,mvaEdb?.z},''comp'');',ic);
hca.YLabel.String = {'JxB_{N}','[mV/m]'};
irf_legend(hca,{'JxB','E'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,mvaJ);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_L','J_M','J_N'},[0.95 0.95]);

if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end

h(1).Title.String = irf_ssub('MMS ?: E dot B = 0',ic);
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

