%% General stuff
c_eval('ExB? = dslE?brst.cross(dmpaB?brst.resample(dslE?brst))/dmpaB?brst.resample(dslE?brst).abs2*1e3;');


%% Overview
c_eval('Pall? = irf.ts_scalar(Pe_psd?.time,[Pe_psd?.data(:,1) Pe_psd?.data(:,4) Pe_psd?.data(:,5) Pe_psd?.data(:,2) Pe_psd?.data(:,6) Pe_psd?.data(:,3)]);')
c_eval('[PeXXp?,~,~,PeYYp?,~,PeZZp?] = mms.rotate_tensor_fac(Pall?,dmpaB?brst);')
c_eval('Pe?perp = (PeXXp?+PeYYp?)/2; Pe?par = PeZZp?;')

c_eval('Tall? = irf.ts_scalar(Te_psd?.time,[Te_psd?.data(:,1) Te_psd?.data(:,4) Te_psd?.data(:,5) Te_psd?.data(:,2) Te_psd?.data(:,6) Te_psd?.data(:,3)]);')
c_eval('[TeXXp?,~,~,TeYYp?,~,TeZZp?] = mms.rotate_tensor_fac(Tall?,dmpaB?brst);')
c_eval('Te?perp = (TeXXp?+TeYYp?)/2; Te?par = TeZZp?;')

%%

tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:05.00Z'); % magnetosphere-magnetosheath-magnetosphere
h = irf_plot(11);

ic  =4;
% Magnetic field
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.7]);



hca = irf_panel('brst n');
set(hca,'ColorOrder',mms_colors('xyza'))
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
%c_eval('irf_plot(hca,{ne_psd1.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
c_eval('irf_plot(hca,{ne_psd1.tlim(tint)});',ic); 
%hca.YLabel.String = {'n','(cm^{-3})'};
ylabel(hca,{'n','(cm^{-3})'},'interpreter','tex')
%hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('xyza'))
%irf_legend(hca,{'n_e','n_i'},[1.01 0.5]);

hca = irf_panel('brst vi');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{vi?_lowres.tlim(tint).x,vi?_lowres.tlim(tint).y,vi?_lowres.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'v_i','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[1.01 0.7]);

hca = irf_panel('brst ve');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{ve?brst.tlim(tint).x,ve?brst.tlim(tint).y,ve?brst.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'v_e','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[1.01 0.7]);

hca = irf_panel('j curl');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{jbrst.x,jbrst.y,jbrst.z},'comp');
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'J_x','J_y','J_z'},[1.01 0.7]);

hca=irf_panel('idist');
c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{i,OMNI}','(eV)'},'Interpreter','tex');

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{e,OMNI}','(eV)'},'Interpreter','tex');

clim = [-1 1];
ylim = [10 1e3];

if 1  
  hca = irf_panel(irf_ssub('edistparperp ?',ic));
  c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
  hold(hca,'on')  
  hold(hca,'off')    
  set(hca,'yscale','log');
  hca.CLim = clim;
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hca.YLim = ylim;
  ylabel(hca,{'E_e','(eV)'});
  colormap(hca,cn.cmap('bluered3'))

  if 1 % Colorbar labels
    hcb.YTick = 0.6*[-1 1];
    hcb.YTickLabel = {'f_{perp}','f_{a/par}'};
    %hcb.YLabel.String = irf_ssub('mms {?}',ic);
    hcb.YLabel.String = ' ';
    %hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
  end
end
if 1
  hca = irf_panel(irf_ssub('edistparapar ?',ic));
  c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)  
  set(hca,'yscale','log');
  hca.CLim = clim;
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %hca.YLim = ylim;
  ylabel(hca,{'E_e','(eV)'});
  colormap(hca,cn.cmap('bluered3'))

  if 1 % Colorbar labels
    hcb.YTick = 0.6*[-1 1];
    hcb.YTickLabel = {'f_{apar}','f_{par}'};
    %hcb.YLabel.String = irf_ssub('mms {?}',ic);
    hcb.YLabel.String = ' ';
    %hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
  end
end

hca = irf_panel(irf_ssub('T?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))  
%c_eval('irf_plot(hca,{Teperp?_lowres.tlim(tint),Tepar?_lowres.tlim(tint)},''comp'');',ic)
c_eval('irf_plot(hca,{Te?perp.tlim(tint),Te?par.tlim(tint)},''comp'');',ic)
hca.YLabel.String = {irf_ssub('T_e',ic),'(eV)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'T_{\perp}','T_{||}'},[1.01 0.7]);
  

if 0
hca = irf_panel('j fac');
set(hca,'ColorOrder',mms_colors('xyz'))
c_eval('irf_plot(hca,{jbrst.dot(dmpaB?brst.resample(jbrst)),jtot?.dot(dmpaB?brst.resample(jtot?))},''comp'');',ic)
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'J_{curl}','J_{moments}'},[1.01 0.7]);
end
% Electric field
hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[1.01 0.7]);



irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'))
irf_zoom(h([1:5 10:11]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end

%% Momentum plot
units = irf_units;
c_eval('vexB?mVm = vexB?*1e-3; vexB?.units = ''mV/m'';')
%c_eval('vixB?mVm = vixB?*1e-3; vixB?.units = ''mV/m'';')
gradPene = gradP/avne/units.e*1e-6*1e-9; gradPene.units = 'mV/m';
gradPene = gradPene*1;
Efactor = 0;
%c_eval('jxBne? = jxB/avNe/e*1e-21; jxBne?.units = ''mV/m'';')
h = irf_plot(9);

ic  = 4;
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('x terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.x-Efactor,1*vexB?mVm.x,1*gradPene.x},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'X-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'v_{e}xB','\nabla \cdot P_e/ne'},[0.95 0.95]);

hca = irf_panel('x terms 2');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.x-Efactor,-1*vexB?mVm.x.resample(dslE?brst)-1*gradPene.x.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'X-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);

hca = irf_panel('y terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.y-Efactor,1*vexB?mVm.y.resample(dslE?brst),1*gradPene.y.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'Y-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'v_{e}xB','\nabla \cdot P_e/ne'},[0.95 0.95]);

hca = irf_panel('y terms 2');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.y-Efactor,-1*vexB?mVm.y.resample(dslE?brst)-1*gradPene.y.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'Y-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);

hca = irf_panel('z terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.z-Efactor,1*vexB?mVm.z.resample(dslE?brst),1*gradPene.z.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'Z-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'v_{e}xB','\nabla \cdot P_e/ne'},[0.95 0.95]);

hca = irf_panel('z terms 2');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.z-Efactor,-1*vexB?mVm.z.resample(dslE?brst)-1*gradPene.z.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'Z-terms','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E-?',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);


if 0
hca = irf_panel(irf_ssub('E-?+ vxB',Efactor));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst),dslE?brst.y-Efactor+vexB?mVm.y.resample(dslE?brst),dslE?brst.z-Efactor+vexB?mVm.z.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time)
hca.YLabel.String = {irf_ssub('E-?+v_exB',Efactor),'(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'x','y','z'},[0.95 0.95]);
end
if 0
  hca = irf_panel(irf_ssub('P?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{PeXXp?.tlim(tint),PeZZp?.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('P',ic),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{\perp}','P_{||}'},[0.95 0.95]);
end
if 1
  hca = irf_panel(irf_ssub('T?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{TeXXp?.tlim(tint),TeZZp?.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('T',ic),'(eV)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'T_{\perp}','T_{||}'},[0.95 0.95]);
end

hca = irf_panel('brst n');
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{ne1_lowres.tlim(tint),ne2_lowres.tlim(tint),ne3_lowres.tlim(tint),ne4_lowres.tlim(tint)},'comp');
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

irf_zoom(h,'x',irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:32.00Z'))
irf_zoom(h,'y')
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Diffusion region
h = irf_plot(8);

ic  = 4;

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brstRemOff.z.tlim(tint),dmpaB2brstRemOff.z.tlim(tint),dmpaB3brstRemOff.z.tlim(tint),dmpaB4brstRemOff.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','(nT)'};
%set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);
set(hca,'ColorOrder',mms_colors('1')); irf_legend(hca,{'mms 1'},[1.01 0.99]);
set(hca,'ColorOrder',mms_colors('2')); irf_legend(hca,{'mms 2'},[1.01 0.79]);
set(hca,'ColorOrder',mms_colors('3')); irf_legend(hca,{'mms 3'},[1.01 0.59]);
set(hca,'ColorOrder',mms_colors('4')); irf_legend(hca,{'mms 4'},[1.01 0.39]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brstRemOff.y.tlim(tint),dmpaB2brstRemOff.y.tlim(tint),dmpaB3brstRemOff.y.tlim(tint),dmpaB4brstRemOff.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);

if 0 % scPot, 4sc
  hca = irf_panel('scPot 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-1*P1brst.tlim(tint),-1*P2brst.tlim(tint),-1*P3brst.tlim(tint),-1*P4brst.tlim(tint)},'comp');
  hca.YLabel.String = {'sc Pot','(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);
end


hca = irf_panel('brst E 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dslE1brst.tlim(tint).x,dslE2brst.tlim(tint).x,dslE3brst.tlim(tint).x,dslE4brst.tlim(tint).x},'comp');
hca.YLabel.String = {'E_x','(mV/m)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);


hca = irf_panel('brst ve 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{ve1brst.tlim(tint).z,ve2brst.tlim(tint).z,ve3brst.tlim(tint).z,ve4brst.tlim(tint).z},'comp');
hca.YLabel.String = {'v_e_z','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);


hca = irf_panel('x terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('lines = irf_plot(hca,{dslE?brst.x-Efactor,1*vexB?mVm.x,1*gradPene.x,-1*vexB?mVm.x.resample(dslE?brst)-1*gradPene.x.resample(dslE?brst)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'E_x','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E',Efactor),'v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
%lines.Children(6).LineWidth = 1.5;
%lines.Children(6).LineStyle = '--';
hca.YLim = [-10 10];

if 0
  hca = irf_panel('x terms 2');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.x-Efactor,-1*vexB?mVm.x.resample(dslE?brst)-1*gradPene.x.resample(dslE?brst)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_x','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{irf_ssub('E',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
  irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
end

hca = irf_panel('y terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('lines = irf_plot(hca,{dslE?brst.y-Efactor,1*vexB?mVm.y.resample(dslE?brst),1*gradPene.y.resample(dslE?brst),-1*vexB?mVm.y.resample(dslE?brst)-1*gradPene.y.resample(dslE?brst)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'E_y','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E',Efactor),'v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
%lines.Children(6).LineWidth = 1.5;
hca.YLim = [-10 10];

if 0
  hca = irf_panel('y terms 2');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.y-Efactor,-1*vexB?mVm.y.resample(dslE?brst)-1*gradPene.y.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{irf_ssub('E',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
  irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
end
% E
if 1 % wavelet
  hca=irf_panel('E wavelet');
  c_eval('irf_spectrogram(hca,wavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
else % pfft
  hca=irf_panel('E pfft');
  c_eval('irf_spectrogram(hca,pfftE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end
% B
if 1 % wavelet
  hca=irf_panel('B wavelet');
  c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
else % pfft
  hca=irf_panel('B pfft');
  c_eval('irf_spectrogram(hca,pfftB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end

irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h(1:6),'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end

%% Diffusion region in LMN coordinates
units = irf_units;
% XX, YY, ZZ, XY, XZ, YZ. 
% Current from 4sc magnetic field: assuming GSE and DMPA are the same coordinate system.
c_eval('gseR?brsttime = gseR?.resample(dmpaB?brstRemOff);',1:4)
[jbrst,divBbrst,Bbrst,jxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','dmpaB?brstRemOff');
jbrst = irf.ts_vec_xyz(jbrst.time,jbrst.data);
jbrst.data = jbrst.data*1e9; jbrst.units = 'nAm^{-2}';
jbrst.time = EpochTT(jbrst.time); jbrst.name = '4sc current density';
%%
if 0 % Calculate pressure gradient
  %% with recalculated moments
  calPf1=1;calPf2=1;calPf3=1;calPf4=1;
  c_eval(['P? = {irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,1))*calPf?,',... % xx
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,4))*calPf?,',... % xy
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,5))*calPf?;',... % xz
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,4))*calPf?,',... % yx
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,2))*calPf?,',... % yy
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,6))*calPf?;',... % yz
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,5))*calPf?,',... % zx
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,6))*calPf?,',... % zy
                'irf.ts_scalar(Pe_psd?.time,Pe_psd?.data(:,3))*calPf?};']);
  gradP = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,P1,P2,P3,P4);
  avne = (ne_psd1.resample(ne_psd1.time) + ne_psd2.resample(ne_psd1.time) + ne_psd3.resample(ne_psd1.time) + ne_psd4.resample(ne_psd1.time))/4;
  units = irf_units;
  gradPene = gradP/avne/units.e*1e-6*1e-9; gradPene.units = 'mV/m'; gradPene.name = 'grad P/ne'
else
  %% old way
  calPf1 = 1; calPf2 = 1; calPf3 = 1.08; calPf4 = 1;
  c_eval('P? = {Pexx?brst*calPf?,Pexy?brst*calPf?,Pexz?brst*calPf?;Pexy?brst*calPf?,Peyy?brst*calPf?,Peyz?brst;Pexz?brst*calPf?,Peyz?brst*calPf?,Pezz?brst*calPf?};')

  gradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,P1,P2,P3,P4);    
  avne = (ne1brst*calPf1 + ne2brst.resample(ne1brst.time)*calPf2 + ne3brst.resample(ne1brst.time)*calPf3 + ne4brst.resample(ne1brst.time)*calPf4)/4;
  gradPene = gradPe/avne/units.e*1e-6*1e-9; gradPene_b.units = 'mV/m'; gradPene.name = 'grad P/ne'
end

c_eval('vexB? = Ve_psd?.resample(dmpaB?brst).cross(dmpaB?brst);')
c_eval('vexB?mVm = vexB?*1e-3; vexB?.units = ''mV/m'';')

avB = (dmpaB1brstRemOff + dmpaB2brstRemOff.resample(dmpaB1brst) + dmpaB3brstRemOff.resample(dmpaB1brst) + dmpaB4brstRemOff.resample(dmpaB1brst))/4;
jxB = jbrst.resample(avB).cross(avB);
jxBmVm = jxB*1e-9*1e-9/(avne*1e6)/units.e*1e3; jxBmVm.units = 'mV/m'; 

%c_eval('vexB? = vi?brst.resample(dmpaB?brst).cross(dmpaB?brst);')
%c_eval('vexB?mVm = vexB?*1e-3; vexB?.units = ''mV/m'';')

%%
tint = irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z');
[out,l,v] = irf_minvar(dmpaB1brst.tlim(tint));

L = v(1,:);
M = v(2,:);
N = v(3,:);
%%
c_eval('plB? = irf.ts_vec_xyz(dmpaB?brstRemOff.time,[dmpaB?brstRemOff.dot(L).data dmpaB?brstRemOff.dot(M).data dmpaB?brstRemOff.dot(N).data]);')
c_eval('plE? = irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(L).data dslE?brst.dot(M).data dslE?brst.dot(N).data]);')
c_eval('plve? = irf.ts_vec_xyz(Ve_psd?.time,[Ve_psd?.dot(L).data Ve_psd?.dot(M).data Ve_psd?.dot(N).data]);')
c_eval('vi?brst.name = ''vi brst'';plvi? = irf.ts_vec_xyz(vi?brst.time,[vi?brst.dot(L).data vi?brst.dot(M).data vi?brst.dot(N).data]);')
plGradPene = irf.ts_vec_xyz(gradPene.time,[gradPene.dot(L).data gradPene.dot(M).data gradPene.dot(N).data]);
pljxB = irf.ts_vec_xyz(jxBmVm.time,[jxBmVm.dot(L).data jxBmVm.dot(M).data jxBmVm.dot(N).data]);
c_eval('plvexB?mVm = irf.ts_vec_xyz(vexB?mVm.time,[vexB?mVm.dot(L).data vexB?mVm.dot(M).data vexB?mVm.dot(N).data]);')
avE = (plE1+ plE2.resample(plE1) + plE3.resample(plE1) + plE4.resample(plE1))/4;
avvexB = (plvexB1mVm + plvexB2mVm.resample(plvexB1mVm) + plvexB3mVm.resample(plvexB1mVm) + plvexB4mVm.resample(plvexB1mVm))/4;
c_eval('plvixB?mVm = plvi?.resample(plB?).cross(plB?)*1e-3;')
avvixB = (plvixB1mVm + plvixB2mVm.resample(plvixB1mVm) + plvixB3mVm.resample(plvixB1mVm) + plvixB4mVm.resample(plvixB1mVm))/4;

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
%tint = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(7);

ic  = 1;

hca = irf_panel('BL');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plB1.x.tlim(tint),plB2.x.tlim(tint),plB3.x.tlim(tint),plB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
%set(hca,'ColorOrder',mms_colors('1')); irf_legend(hca,{'mms 1'},[0.98 0.99]);
%set(hca,'ColorOrder',mms_colors('2')); irf_legend(hca,{'mms 2'},[0.98 0.79]);
%set(hca,'ColorOrder',mms_colors('3')); irf_legend(hca,{'mms 3'},[0.98 0.59]);
%set(hca,'ColorOrder',mms_colors('4')); irf_legend(hca,{'mms 4'},[0.98 0.25]);

hca = irf_panel('BM');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plB1.y.tlim(tint),plB2.y.tlim(tint),plB3.y.tlim(tint),plB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{M}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);

if 1
hca = irf_panel('BN');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plB1.z.tlim(tint),plB2.z.tlim(tint),plB3.z.tlim(tint),plB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{N}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);
end

if 0 % scPot, 4sc
  hca = irf_panel('scPot 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-1*P1brst.tlim(tint),-1*P2brst.tlim(tint),-1*P3brst.tlim(tint),-1*P4brst.tlim(tint)},'comp');
  hca.YLabel.String = {'sc Pot','(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))
%  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);
end

 
if 0
  hca = irf_panel('brst EM 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{plE1.tlim(tint).y,plE2.tlim(tint).y,plE3.tlim(tint).y,plE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
%  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);
end

hca = irf_panel('brst ve L 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plve1.tlim(tint).x,plve2.tlim(tint).x,plve3.tlim(tint).x,plve4.tlim(tint).x},'comp');
hca.YLabel.String = {'v_e_L','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);

hca = irf_panel('brst EN 4sc');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plE1.tlim(tint).z,plE2.tlim(tint).z,plE3.tlim(tint).z,plE4.tlim(tint).z},'comp');
hca.YLabel.String = {'E_N','(mV/m)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[1.01 0.9]);


if 0
  hca = irf_panel('N terms');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('lines = irf_plot(hca,{plE?.z,1*plvexB?mVm.z,1*plGradPene.z,-1*plvexB?mVm.z.resample(plE?)-1*plGradPene.z.resample(plE?)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.05]);
  irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1
  %%
   hca = irf_panel('N terms');
  set(hca,'ColorOrder',mms_colors('xyzb'))
  %irf_plot(hca,{avE.z,avvexB.z,plGradPene.z,-1*avvexB.z.resample(avE)-1*plGradPene.z.resample(avE)},'comp'); 
  irf_plot(hca,{avE.z,-1*avvexB.z,-1*plGradPene.z,-1*avvexB.z.resample(avE)-1*plGradPene.z.resample(avE)},'comp'); 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzb'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end

if 0
hca = irf_panel('M terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('lines = irf_plot(hca,{plE?.y,1*plvexB?mVm.y,1*plGradPene.y,-1*plvexB?mVm.y.resample(plE?)-1*plGradPene.y.resample(plE?)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'E_M','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95],'fontsize',12);
irf_legend(hca,irf_ssub('4 sc av',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
%lines.Children(6).LineWidth = 1.5;
%lines.Children(6).LineStyle = '--';
hca.YLim = [-10 10];
end
if 0
  hca = irf_panel('av E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{avE.x,avE.y,avE.z},'comp') 
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.95 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);
  hca.YLim = [-10 10];

  hca = irf_panel('jxB terms');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{pljxB.x,pljxB.y,pljxB.z},'comp')
  hca.YLabel.String = {'jxB/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.95 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('curlometer',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);
  hca.YLim = [-10 10];


  hca = irf_panel('vixB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{-1*avvixB.x,-1*avvixB.y,-1*avvixB.z},'comp') 
  hca.YLabel.String = {'-v_ixB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.95 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);
  hca.YLim = [-10 10];

  hca = irf_panel('av gradPe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{-1*plGradPene.x,-1*plGradPene.y,-1*plGradPene.z},'comp') 
  hca.YLabel.String = {'-\nabla\cdot P_e/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.95 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);
  hca.YLim = [-10 10];
else
  hca = irf_panel('ohms law N');
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_plot(hca,{avE.z,pljxB.z,-1*avvixB.z,-1*plGradPene.z.resample(pljxB),pljxB.z-1*avvixB.z.resample(pljxB)-1*plGradPene.z.resample(pljxB)},'comp') 
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color',mms_colors(irf_ssub('?',ic)),'fontsize',12);

  hca.YLim = [-10 10];

end


 

irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h,'x',irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z'))
irf_zoom(h,'x',irf.tint('2015-10-16T10:33:23.00Z/2015-10-16T10:33:34.00Z'))
irf_zoom(h,'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%add_length_on_top(h(1),20,0.5)
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.02 0.98],'color',[0 0 0])
  %h(ii).Position(3) = h(ii).Position(3)*0.95;
end

%%
%%
if 0
  hca = irf_panel('x terms 2');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.x,-1*vexB?mVm.x.resample(dslE?brst)-1*gradPene.x.resample(dslE?brst)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_x','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
  irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
end

hca = irf_panel('y terms');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('lines = irf_plot(hca,{dslE?brst.y,1*vexB?mVm.y.resample(dslE?brst),1*gradPene.y.resample(dslE?brst),-1*vexB?mVm.y.resample(dslE?brst)-1*gradPene.y.resample(dslE?brst)},''comp'');',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
% ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
hca.YLabel.String = {'E_y','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{irf_ssub('E',Efactor),'v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
%lines.Children(6).LineWidth = 1.5;
hca.YLim = [-10 10];

if 0
  hca = irf_panel('y terms 2');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.y-Efactor,-1*vexB?mVm.y.resample(dslE?brst)-1*gradPene.y.resample(dslE?brst)},''comp'')',ic) % ,dslE?.x.resample(vexB?mVm.time)+vexB?mVm.x.resample(vexB?mVm.time)+gradPene?.x.resample(vexB?mVm.time) 
  % ,dslE?brst.x-Efactor+vexB?mVm.x.resample(dslE?brst)+gradPene.x.resample(dslE?brst)
  hca.YLabel.String = {'E_y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{irf_ssub('E',Efactor),'-v_{e}xB-\nabla \cdot P_e/ne'},[0.95 0.95]);
  irf_legend(hca,irf_ssub('mms ?',ic),[1.01 0.95],'color',mms_colors(irf_ssub('?',ic)));
end
% E
if 0 % wavelet
  hca=irf_panel('E wavelet');
  c_eval('irf_spectrogram(hca,wavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('E pfft');
  c_eval('irf_spectrogram(hca,pfftE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end
% B
if 0 % wavelet
  hca=irf_panel('B wavelet');
  c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('B pfft');
  c_eval('irf_spectrogram(hca,pfftB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end

irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h(1:6),'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end

%% Waves

tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:05.00Z'); % magnetosphere-magnetosheath-magnetosphere

%c_eval('plE? = irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(L).data dslE?brst.dot(M).data dslE?brst.dot(N).data]);')
c_eval('plE? = dslE?brst;')

h = irf_plot(8);

ic  =4;
% Magnetic field
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.7]);

% Electric field
hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[1.01 0.7]);

% E
hca=irf_panel('E wavelet'); irf_plot(hca,dslE1brst);
if 0 % wavelet  
  c_eval('irf_spectrogram(hca,wavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('E pfft');
  c_eval('irf_spectrogram(hca,pfftE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end
% B
  hca=irf_panel('B wavelet'); irf_plot(hca,dmpaB1scm);
if 0 % wavelet
  c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('B pfft');
  c_eval('irf_spectrogram(hca,pfftB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end

irf_zoom(h(1:4),'x',tint)
%irf_pl_mark(h(1:2),tint)
hca = irf_panel('delete for space'); %delete(hca)    
grid(hca,'off')    
irf_timeaxis(h(4))
    
if 1 %plot lh
  hca = irf_panel('lh phi')
  tint_lh = irf.tint('2015-10-16T10:33:24.00Z',1.5); 
  E = dslE3brst.tlim(tint_lh);
  acB = dmpaB3scm.tlim(tint_lh);
  dcB = dmpaB3.tlim(tint_lh);
  c_eval('n_loc = mean(ne_psd?.tlim(tint_lh).data);',sc)
  c_eval('Te_loc = mean(Te?perp.tlim(tint_lh).data);',sc)

  angles = 0:2:359;
  f_highpass = 0.5;

  % Propagation direction
  [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir({dcB,acB},E,angles,f_highpass);
  i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
  direction=x(i_dir,:);

  % Velocity and density
  mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3.mean

  v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
  phiB_scale = B0*1e-18/mu0/e/n; % km/s
  v = v_approx*linspace(0.1,2.5,30);

  [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
  i_v=find(corr_v==min(corr_v));
  velocity=v(i_v);
  tsBz = irf.ts_scalar(irf_time(Bz(:,1),'epoch>epochunix'),Bz(:,2));
  tsIntEdt = irf.ts_scalar(irf_time(intEdt(:,1),'epoch>epochunix'),intEdt(:,1+i_dir));
  irf_plot(hca,{tsBz*phiB_scale,tsIntEdt*velocity},'comp')
  %irf_legend(hca,{'B_z*scale','\int Edt*v'},[0.98 0.95])
  irf_legend(hca,{'\phi_{B}','\phi_{E}'},[0.98 0.95])
  hca.YLabel.String = '\Phi [V]'; 
  irf_legend(hca,{['\omega > ' num2str(flim,'%.2f')]},[0.02 0.1])
  irf_zoom(hca,'x',tint_lh)
  %irf_zoom(hca,'y',[-300 300])    
  axtop = add_length_on_top(hca,velocity,0.5);
  axtop.XLabel.String = 'Length [km]';    
  % add B labels on right
  yratio = phi_B(100,2)/Bz(100,2);            
  axtop.YLim = hca.YLim/yratio;
  axtop.YLabel.String = '\delta B_{||} [nT]'; 
  axtop.YTickLabelMode = 'auto';
  axtop.YTickMode = 'auto';
  % mark different physical lengths
  irf_plot_axis_align
  irf_plot_zoomin_lines_between_panels(h(4),h(6))
  hca = irf_panel('delete for space'); delete(hca) 
    
  % Compare Ek and Bz amplitudes
  maxEk = max(dEk(:,2));
  maxEn = max(dEn(:,2));
  maxBz = max(Bz(:,2));
  disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
end   

hca = irf_panel('delete for space 2'); %delete(hca)    
grid(hca,'off')    
irf_timeaxis(h(6))

c_eval('Epar = irf.ts_scalar(dslE?brst.time,dslE?brst.dot(dmpaB1brst.resample(dslE?brst)).data)',ic)



irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h(1:6),'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end

%% Lower hybrid waves

%% Electron outflow
  %%
  ic =1;
  %tint = irf.tint('2015-10-16T10:33:29.40Z',1);
  for ii = 11%:17;
  tint = irf.tint('2015-10-16T10:33:30.00Z',0.01); tint = tint(1);  
  tint = tint+0.03*(ii-1);
  tint = tint(1);
  for ii = 1:4
    h2(ii) = subplot(2,2,ii);
  end
  ic = 3;
  %ii= ii+1;
  %if ii>10; break; end
  %tint = tint + 0.1; 
  %tint = tint(ii);
 
  
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [-2 5.2];
  palim = [1e-2 1e6];
  skymapEnergy = [178];
  %skymapEnergy = [138];
  
  c_eval('dist = desDist?;',ic)

  %c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  %hatVi0 = double(irf_norm(Vi0));
  %c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  %hatVe0 = double(irf_norm(Ve0));

  % Get mean magnetic field direction
  %c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  c_eval('B0 = mean(dmpaB?brst.resample(dmpaB?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  c_eval('E0 = mean(dslE?brst.tlim(tint+[-0.01 0.01]).data);',ic); 
  hatB0 = double(irf_norm(B0));   
  hatE0 = double(irf_norm(E0));   
  ExB0 = cross(E0,B0)/norm(B0)/norm(B0);
  hatExB0 = ExB0/norm(ExB0);
  vectors = {hatB0,'B';[0 0 1],'Z';hatExB0,'ExB';hatE0,'E'};%0;hatVe0,'V_e'};
  
   
  z = hatB0;
  x = cross(cross(z,hatExB0),z); x = x/norm(x);  
  y = cross(z,x); y = y/norm(y);
  %x = L; y = M; z = N;
  %mva
  
  isub = 1;
  
  % Plot psd 0 90 180
  hca = h2(isub); isub = isub + 1;
  c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint+0.03*[-1 1],''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
  %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
  
  if 0
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
  end
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors(1,1:2),''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
 % hca.CLim = [0 5000];
  
  % Plot project ion onto a plane
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)

  if 1
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  else
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  end
  %cn.print([irf_ssub('mms?_',ic) irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);

  end

%% LABELS ETC ADDED

  %tint = irf.tint('2015-10-16T10:33:29.40Z',1);
  for ii = 6%:17;
  tint = irf.tint('2015-10-16T10:33:30.30Z',0.01); tint = tint(1);  
  tint = tint+-0.03*(ii-1);
  tint = tint(1);
  for ii = 1:4
    h2(ii) = subplot(2,2,ii);
  end
  ic = 4;
  %ii= ii+1;
  %if ii>10; break; end
  %tint = tint + 0.1; 
  %tint = tint(ii);
 
  
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [-2 5.2];
  palim = [1e-2 1e6];
  skymapEnergy = [178];
  %skymapEnergy = [138];
  
  c_eval('dist = desDist?;',ic)

  %c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  %hatVi0 = double(irf_norm(Vi0));
  %c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  %hatVe0 = double(irf_norm(Ve0));

  % Get mean magnetic field direction
  %c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  c_eval('B0 = mean(dmpaB?brst.resample(dmpaB?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  c_eval('E0 = mean(dslE?brst.tlim(tint+[-0.01 0.01]).data);',ic); 
  hatB0 = double(irf_norm(B0));   
  hatE0 = double(irf_norm(E0));   
  ExB0 = cross(E0,B0)/norm(B0)/norm(B0);
  hatExB0 = ExB0/norm(ExB0);
  vectors = {hatB0,'B';[0 0 1],'Z';hatExB0,'ExB';hatE0,'E'};%0;hatVe0,'V_e'};
  
   
  z = hatB0;
  x = cross(cross(z,hatExB0),z); x = x/norm(x);  
  y = cross(z,x); y = y/norm(y);
  %x = L; y = M; z = N;
  %mva
  
  
  isub = 1;
  
  % Plot psd 0 90 180
  hca = h2(isub); isub = isub + 1;
  c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint+0.03*[-1 1],''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
  %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
  hca.Position(3)=hca.Position(3)*0.8;
  
  if 0
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
  end
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors(1,1:2),''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
 % hca.CLim = [0 5000];
  hca.Position(3)=hca.Position(3)*1.25;
  hca.Position(1)=hca.Position(1)-0.08;
  
  % Plot project ion onto a plane
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;x;-z],'vlabel',{'E','ExB','B'},'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  %mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;x;-z],'vlabel',{'E','ExB','B'},'elevationlim',elevlim,'vlim',vlim,'vectors',{L,'L';M,'M';N,'N'},'clim',projclim);
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'vlabel',{'v_{ExB dir}','v_{E dir}','v_{B dir}'},'elevationlim',elevlim,'vlim',vlim,'vectors',{L,'L';M,'M';N,'N'},'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)

  if 1
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'vlabel',{'B','ExB','E'},'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  %mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'vlabel',{'B','ExB','E'},'elevationlim',elevlim,'vlim',vlim,'vectors',{L,'L';M,'M';N,'N'},'clim',projclim);
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;z;-y],'vlabel',{'v_{ExB dir}','v_{B dir}','v_{E dir}'},'elevationlim',elevlim,'vlim',vlim,'vectors',{L,'L';M,'M';N,'N'},'clim',projclim);%{L,'L';M,'M';N,'N'},'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  else
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  end
  
  h2(1).Title.String = h2(2).Title.String{1};
  h2(2).Title.String{1} = [];
  %cn.print([irf_ssub('mms?_',ic) irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);

  irf_legend(h2(1),{'(a)'},[0.02 0.98],'color','black')
  irf_legend(h2(2),{'(b)'},[0.02 0.98],'color','white')
  irf_legend(h2(3),{'(c)'},[0.02 0.98],'color','white')
  irf_legend(h2(4),{'(d)'},[0.02 0.98],'color','white')
  end







%% Check knee distributions
tint = irf.tint('2015-10-16T10:33:30.320Z',0.01); tint = tint(1);
vlim = 15*1e3;
elevlim = 20;
correctBin = 0;


c_eval('dist = desDist?;',ic)

%c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
%hatVi0 = double(irf_norm(Vi0));
c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatVe0 = double(irf_norm(Ve0));

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));

% Initialize figure
for ii = 1:8; h(ii) = subplot(2,4,ii); end

isub = 1;

% Plot projection onto a plane
hca = h(isub);  isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; 0 0 1; 0 1 0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 1 0; 0 0 1],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; cross(hatB0,[0 1 0])],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; hatB0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);
