  % Plots single satellite overview of diffusion region using burst data
% First run mms_Oct16.load_brst_data

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z');
tref = tint+[1 0];
ic = 1;

c_eval('Pref? = Pe?_lowres.tlim(tref).data(1,:);')
c_eval('Pe?rel1 = Pe?_lowres+Pref1-Pref?;') % Pa
%irf_plot({Pe2rel1,Pe2_lowres},'comp'); legend('Mod','Orig')
%figure;irf_plot('Pe?rel1','comp')
%legend({'mms1','mms2','mms3','mms4'})
%% Calculate grad P.abs

c_eval('r? = [Pe1rel1.time.epochUnix gseR?.resample(Pe1rel1.time).data Pe?rel1.abs.resample(Pe1rel1.time).data/3];')
for ii = 1:size(r1,1)
  gradPhi(ii,1)=r1(ii,1);
  gradPhi(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
gradP = irf.ts_vec_xyz(irf_time(gradPhi(:,1),'epoch>utc'),gradPhi(:,2:4));
%%

legends4SC = {'1','2','3','4'};
h = irf_plot(11);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('brst E');
irf_plot(hca,{dslE1brst.x,dslE2brst.x,dslE3brst.x,dslE4brst.x},'comp');
hca.YLabel.String = {'E_{X,DSL}','[mV/m]'};
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst n');
irf_plot(hca,{ne1_lowres,ne2_lowres,ne3_lowres,ne4_lowres},'comp');,
hca.YLabel.String = {'n','[cm^{-3}]'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Tx');
irf_plot(hca,{Te1_lowres.x,Te2_lowres.x,Te3_lowres.x,Te4_lowres.x},'comp');
hca.YLabel.String = {'T_{e,x}','(eV)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Ty');
irf_plot(hca,{Te1_lowres.y,Te2_lowres.y,Te3_lowres.y,Te4_lowres.y},'comp');
hca.YLabel.String = {'T_{e,y}','(eV)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Tz');
irf_plot(hca,{Te1_lowres.z,Te2_lowres.z,Te3_lowres.z,Te4_lowres.z},'comp');
hca.YLabel.String = {'T_{e,z}','(eV)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Px');
irf_plot(hca,{Pe1rel1.x*1e9,Pe2rel1.x*1e9,Pe3rel1.x*1e9,Pe4rel1.x*1e9},'comp');
hca.YLabel.String = {'P_{e,x}','(nPa)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Py');
irf_plot(hca,{Pe1rel1.y*1e9,Pe2rel1.y*1e9,Pe3rel1.y*1e9,Pe4rel1.y*1e9},'comp');
hca.YLabel.String = {'P_{e,y}','(nPa)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('brst Pz');
irf_plot(hca,{Pe1rel1.z*1e9,Pe2rel1.z*1e9,Pe3rel1.z*1e9,Pe4rel1.z*1e-9},'comp');
hca.YLabel.String = {'P_{e,z}','(nPa)'};
hca.YScale = 'lin';
irf_legend(hca,legends4SC,[0.95 0.95]);


hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst);',ic);
hca.YLabel.String = {'v_e','(km/s)'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);


irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h,'y')
if 0
  %%
  tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:55.00Z');
  irf_zoom(h,'x',tint)
  irf_zoom(h,'y')
end

%%

legends4SC = {'1','2','3','4'};
h = irf_plot(7);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('brst E');
irf_plot(hca,{dslE1.x,dslE2.x,dslE3.x,dslE4.x},'comp');
hca.YLabel.String = {'E_{X,DSL}','[mV/m]'};
irf_legend(hca,legends4SC,[0.95 0.95]);

hca = irf_panel('JxB');
e = 1.6022e-19;
avNe = (ne1brst + ne2brst.resample(ne1brst.time) + ne3brst.resample(ne1brst.time) + ne4brst.resample(ne1brst.time))/4;
avB = (dmpaB1.resample(dmpaB1brst.time) + dmpaB2.resample(dmpaB1brst.time) + dmpaB3.resample(dmpaB1brst) + dmpaB4.resample(dmpaB1brst))/4;
jxB = j.cross(avB.tlim(tint).resample(j.time)); 
unit_factor = 1e-21;
 % going to mV/m
irf_plot(hca,jxB/avNe/e*unit_factor); % nAm^{-2}*nT/cm^{-3} = 10^-18*10^6 m^3 Am^{-2}*T/e = 10^-12 AmT/e = 10^-12/e V/m = 10^-9/e mV/m   mV/m
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'JxB/ne','[nA m^{-2} T]'},'Interpreter','tex');
irf_legend(hca,{'jxB_X','jxB_y','jxB_z'},[0.95 0.95]);

hca = irf_panel('gradPe');
unit_factor = 1e-6;
irf_plot(hca,gradP/avNe/e*unit_factor); % Pa/km*10^6/m^3/e = 10^6/m^3*Pa/(m 10^3)/e = 10^6*10^-3*Pa/m^4/e = 10^3*Pa/m^4/e
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'grad Pe/ne','[Pa km^{-1}]'},'Interpreter','tex');
irf_legend(hca,{'x','y','z'},[0.95 0.95]);


hca = irf_panel('vexB');
unit_factor = 1e-3;
c_eval('vexB? = ve?brst.cross(dmpaB?brst.resample(ve?brst.time));')
irf_plot(hca,vexB1*unit_factor); % Pa/km*10^6/m^3/e = 10^6/m^3*Pa/(m 10^3)/e = 10^6*10^-3*Pa/m^4/e = 10^3*Pa/m^4/e
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'v_e x B','[]'},'Interpreter','tex');
irf_legend(hca,{'x','y','z'},[0.95 0.95]);

hca = irf_panel('vixB');
unit_factor = 1e-3;
c_eval('vixB? = vi?brst.cross(dmpaB?brst.resample(vi?brst.time));')
irf_plot(hca,vixB1*unit_factor); % Pa/km*10^6/m^3/e = 10^6/m^3*Pa/(m 10^3)/e = 10^6*10^-3*Pa/m^4/e = 10^3*Pa/m^4/e
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'v_i x B','[]'},'Interpreter','tex');
irf_legend(hca,{'x','y','z'},[0.95 0.95]);
%irf_plot({dmpaB1.tlim(tint),j.tlim(tint),j.tlim(tint).cross(dmpaB1.tlim(tint)),gradP*1e3})

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:55.00Z');
  irf_zoom(h,'x',tint)
  irf_zoom(h,'y')
  
%% Electron momentum plot
ic = 4;

c_eval('vexB?mVm = vexB?*1e-3; vexBne?.units = ''mV/m'';')
c_eval('vixB?mVm = vixB?*1e-3; vixBne?.units = ''mV/m'';')
c_eval('gradPene? = gradP/avNe/e*1e-6; gradPene?.units = ''mV/m'';')
c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')

h = irf_plot(8);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('x terms');
c_eval('irf_plot(hca,{dslE?.x,vexB?mVm.x,gradPene?.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'X-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e,1}xB/ne',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of x terms');
c_eval('irf_plot(hca,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time))',ic)
hca.YLabel.String = {'X-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB/ne+grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('y terms');
c_eval('irf_plot(hca,{dslE?.y,vexB?mVm.y,gradPene?.y},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Y-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e,1}xB/ne',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of y terms');
c_eval('irf_plot(hca,dslE?.y.resample(vexB?mVm.time)+vexB?mVm.y.resample(vexB?mVm.time)+gradPene?.y.resample(vexB?mVm.time))',ic)
hca.YLabel.String = {'Y-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB/ne+grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('z terms');
c_eval('irf_plot(hca,{dslE?.z,vexB?mVm.z,gradPene?.z},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e,1}xB/ne',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of z terms');
c_eval('irf_plot(hca,dslE?.z.resample(vexB?mVm.time)+vexB?mVm.z.resample(vexB?mVm.time)+gradPene?.z.resample(vexB?mVm.time))',ic)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB/ne+grad(P_e)/ne'},[0.95 0.95]);

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
% Compare momentum terms


%% Electron momentum plot
ic = 2;

%c_eval('vexB?mVm = vexB?*1e-3; vexBne?.units = ''mV/m'';')
%c_eval('vixB?mVm = vixB?*1e-3; vixBne?.units = ''mV/m'';')
%%c_eval('gradPene? = gradP/avNe/e*1e-6; gradPene?.units = ''mV/m'';')
%c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')

h = irf_plot(8);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('x terms');
c_eval('irf_plot(hca,{dslE?.x,vexB?.x,gradPene.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'X-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

if 0
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene.x.resample(vexB?.time))',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,{dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene.x.resample(vexB?.time),dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end

hca = irf_panel('y terms');
c_eval('irf_plot(hca,{dslE?.y,vexB?.y,gradPene.y},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Y-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);


if 0
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene.y.resample(vexB?.time))',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,{dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene.y.resample(vexB?.time),dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end



hca = irf_panel('z terms');
c_eval('irf_plot(hca,{dslE?.z,vexB?.z,gradPene.z},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of z terms');
c_eval('irf_plot(hca,dslE?.z.resample(vexB?.time)+vexB?.z.resample(vexB?.time)+gradPene.z.resample(vexB?.time))',ic)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
h(1).Title.String = irf_ssub('MMS ?',ic);
irf_plot_axis_align
% Compare momentum terms

%% Electron momentum plot
ic = 4;

%c_eval('vexB?mVm = vexB?*1e-3; vexBne?.units = ''mV/m'';')
%c_eval('vixB?mVm = vixB?*1e-3; vixBne?.units = ''mV/m'';')
%%c_eval('gradPene? = gradP/avNe/e*1e-6; gradPene?.units = ''mV/m'';')
%c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')

h = irf_plot(8);

hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('1234'))
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{j.x,j.y,j.z},'comp');
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','(nA m^{-2})'},'Interpreter','tex');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('x terms');
set(hca,'ColorOrder',mms_colors('xyz'))
c_eval('irf_plot(hca,{dslE?.x,vexB?.x,gradPene?.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'X-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

if 0
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene?.x.resample(vexB?.time))',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,{dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene?.x.resample(vexB?.time),dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end

hca = irf_panel('y terms');
c_eval('irf_plot(hca,{dslE?.y,vexB?.y,gradPene?.y},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Y-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);


if 0
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene?.y.resample(vexB?.time))',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,{dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene?.y.resample(vexB?.time),dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end



hca = irf_panel('z terms');
c_eval('irf_plot(hca,{dslE?.z,vexB?.z,gradPene?.z},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of z terms');
c_eval('irf_plot(hca,dslE?.z.resample(vexB?.time)+vexB?.z.resample(vexB?.time)+gradPene?.z.resample(vexB?.time))',ic)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
h(1).Title.String = irf_ssub('MMS ?',ic);
irf_plot_axis_align
% Compare momentum terms

%% Generalized Ohm's law
ic = 2;

%c_eval('vexB?mVm = vexB?*1e-3; vexBne?.units = ''mV/m'';')
%c_eval('vixB?mVm = vixB?*1e-3; vixBne?.units = ''mV/m'';')
%%c_eval('gradPene? = gradP/avNe/e*1e-6; gradPene?.units = ''mV/m'';')
%c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')

h = irf_plot(8);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'j_X','j_y','j_z'},[0.95 0.95]);

hca = irf_panel('x terms');
c_eval('irf_plot(hca,{dslE?.x,vexB?.x,gradPene?.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'X-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

if 0
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene?.x.resample(vexB?.time))',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of x terms');
  c_eval('irf_plot(hca,{dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)+gradPene?.x.resample(vexB?.time),dslE?.x.resample(vexB?.time)+vexB?.x.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'X-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end

hca = irf_panel('y terms');
c_eval('irf_plot(hca,{dslE?.y,vexB?.y,gradPene?.y},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Y-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);


if 0
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene?.y.resample(vexB?.time))',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);
else
  hca = irf_panel('sum of y terms');
  c_eval('irf_plot(hca,{dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)+gradPene?.y.resample(vexB?.time),dslE?.y.resample(vexB?.time)+vexB?.y.resample(vexB?.time)},''comp'')',ic)
  hca.YLabel.String = {'Y-terms','(mV/m)'};
  irf_legend(hca,{'E+v_exB+grad(P_e)/ne','E+v_exB'},[0.95 0.95]);
end



hca = irf_panel('z terms');
c_eval('irf_plot(hca,{dslE?.z,vexB?.z,gradPene?.z},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{irf_ssub('E?',ic),irf_ssub('v_{e}xB',ic),'grad(P_e)/ne'},[0.95 0.95]);

hca = irf_panel('sum of z terms');
c_eval('irf_plot(hca,dslE?.z.resample(vexB?.time)+vexB?.z.resample(vexB?.time)+gradPene?.z.resample(vexB?.time))',ic)
hca.YLabel.String = {'Z-terms','(mV/m)'};
irf_legend(hca,{'E+v_exB+grad(P_e)/ne'},[0.95 0.95]);

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
h(1).Title.String = irf_ssub('MMS ?',ic);
irf_plot_axis_align
% Compare momentum terms

