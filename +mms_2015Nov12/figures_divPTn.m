%% Pressure and temperature divergences
% Calibrate density data so that they are at the same level
tref = irf_time('2015-11-12T07:19:20.00Z','utc>EpochTT');
c_eval('ne? = ne?-(ne?.resample(tref).data-ne1.resample(tref).data);',1:4)

gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km';  gseGradPe.name = 'div Pe';
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';  gseGradPi.name = 'div Pi';
gseGradTe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';  gseGradTe.name = 'div Te';
gseGradTi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';  gseGradTi.name = 'div Ti';
gseGradNe = c_4_grad('gseR?','ne?','grad');  gseGradNe.units = 'cc/km';  gseGradNe.name = 'grad Ne';
gseGradNi = c_4_grad('gseR?','ni?','grad');  gseGradNi.units = 'cc/km';  gseGradNi.name = 'grad Ni';

avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
avTe = (gseTe1.trace/3+gseTe2.trace.resample(gseTe1.time)/3+gseTe3.trace.resample(gseTe1.time)/3+gseTe4.trace.resample(gseTe1.time)/3)/4;  avTe.name = '<Te>';
avPe = (gsePe1.trace/3+gsePe2.trace.resample(gsePe1.time)/3+gsePe3.trace.resample(gsePe1.time)/3+gsePe4.trace.resample(gsePe1.time)/3)/4;  avPe.name = '<Pe>';
avNi = (ni1+ni2.resample(ni1.time)+ni3.resample(ni1.time)+ni4.resample(ni1.time))/4; avNi.name = '<ni>';
avTi = (gseTi1.trace/3+gseTi2.trace.resample(gseTi1.time)/3+gseTi3.trace.resample(gseTi1.time)/3+gseTi4.trace.resample(gseTi1.time)/3)/4;  avTi.name = '<Ti>';
avPi = (gsePi1.trace/3+gsePi2.trace.resample(gsePi1.time)/3+gsePi3.trace.resample(gsePi1.time)/3+gsePi4.trace.resample(gsePi1.time)/3)/4;  avPi.name = '<Pi>';

epsNe = gseGradNe/avNe; epsNe.name = 'grad(ne)/ne'; epsNe.units = '1/km';
epsTe = gseGradTe/avTe; epsTe.name = 'grad(Te)/Te'; epsTe.units = '1/km';
epsPe = gseGradPe/avPe; epsPe.name = 'grad(Pe)/Pe'; epsPe.units = '1/km';
epsNi = gseGradNi/avNi; epsNi.name = 'grad(ni)/ni'; epsNi.units = '1/km';
epsTi = gseGradTi/avTi; epsTi.name = 'grad(Ti)/Ti'; epsTi.units = '1/km';
epsPi = gseGradPi/avPi; epsPi.name = 'grad(Pi)/Pi'; epsPi.units = '1/km';

mvaEpsNe = epsNe*lmn'; mvaEpsNe.name = 'grad(ne)/ne'; mvaEpsNe.units = '1/km';
mvaEpsTe = epsTe*lmn'; mvaEpsTe.name = 'grad(Te)/Te'; mvaEpsTe.units = '1/km';
mvaEpsPe = epsPe*lmn'; mvaEpsPe.name = 'grad(Pe)/Pe'; mvaEpsPe.units = '1/km';
mvaEpsNi = epsNi*lmn'; mvaEpsNi.name = 'grad(ni)/ni'; mvaEpsNi.units = '1/km';
mvaEpsTi = epsTi*lmn'; mvaEpsTi.name = 'grad(Ti)/Ti'; mvaEpsTi.units = '1/km';
mvaEpsPi = epsPi*lmn'; mvaEpsPi.name = 'grad(Pi)/Pi'; mvaEpsPi.units = '1/km';


EfromdivPi = avTi.resample(epsNe)*units.eV/units.e*epsNe;
EfromdivPe = avTe*units.eV/units.e*epsNe.resample(avTe);

mvaE_TiepsPi = avTi*units.eV/units.e*mvaEpsPi.resample(avTi);
mvaE_TiepsNe = avTi*units.eV/units.e*mvaEpsNe.resample(avTi);
mvaE_TiepsNi = avTi*units.eV/units.e*mvaEpsNi.resample(avTi);
mvaE_TeepsNe = avTe*units.eV/units.e*mvaEpsNe.resample(avTe);

%% Gradient from 1 sc, using v_cs
CS_normal_velocity = 1*70;
c_eval('ne?_ = ne?.filt(0,10,[],3); dtne? = ne?_.time(2)-ne?_.time(1); dxne? = dtne?*CS_normal_velocity; dne? = ne?_.data(2:end)-ne?_.data(1:end-1); avne? = 0.5*(ne?_.data(2:end)+ne?_.data(1:end-1)); epsNe? = irf.ts_scalar(ne?_.time(1:end-1)+0.5*dtne,dne?./avne?/dxne?); epsNe?.units = ''1/km'';',1:4)
c_eval('Te?_ = mvaTe?.trace/3;      dtTe? = Te?_.time(2)-Te?_.time(1); dxTe? = dtTe?*CS_normal_velocity; dTe? = Te?_.data(2:end)-Te?_.data(1:end-1); avTe? = 0.5*(Te?_.data(2:end)+Te?_.data(1:end-1)); epsTe? = irf.ts_scalar(Te?_.time(1:end-1)+0.5*dtTe,dTe?./avTe?/dxTe?); epsTe?.units = ''1/km'';',1:4)
c_eval('Pe?_ = mvaPe?.trace/3;      dtPe? = Pe?_.time(2)-Pe?_.time(1); dxPe? = dtPe?*CS_normal_velocity; dPe? = Pe?_.data(2:end)-Pe?_.data(1:end-1); avPe? = 0.5*(Pe?_.data(2:end)+Pe?_.data(1:end-1)); epsPe? = irf.ts_scalar(Pe?_.time(1:end-1)+0.5*dtPe,dPe?./avPe?/dxPe?); epsPe?.units = ''1/km'';',1:4)
c_eval('ni?_ = ni?.filt(0,3,[],3);  dtni? = ni?_.time(2)-ni?_.time(1); dxni? = dtni?*CS_normal_velocity; dni? = ni?_.data(2:end)-ni?_.data(1:end-1); avni? = 0.5*(ni?_.data(2:end)+ni?_.data(1:end-1)); epsNi? = irf.ts_scalar(ni?_.time(1:end-1)+0.5*dtni,dni?./avni?/dxni?); epsNi?.units = ''1/km'';',1:4)
c_eval('Ti?_ = mvaTi?.trace/3;      dtTi? = Ti?_.time(2)-Ti?_.time(1); dxTi? = dtTi?*CS_normal_velocity; dTi? = Ti?_.data(2:end)-Ti?_.data(1:end-1); avTi? = 0.5*(Ti?_.data(2:end)+Ti?_.data(1:end-1)); epsTi? = irf.ts_scalar(Ti?_.time(1:end-1)+0.5*dtTi,dTi?./avTi?/dxTi?); epsTi?.units = ''1/km'';',1:4)
c_eval('Pi?_ = mvaPi?.trace/3;      dtPi? = Pi?_.time(2)-Pi?_.time(1); dxPi? = dtPi?*CS_normal_velocity; dPi? = Pi?_.data(2:end)-Pi?_.data(1:end-1); avPi? = 0.5*(Pi?_.data(2:end)+Pi?_.data(1:end-1)); epsPi? = irf.ts_scalar(Pi?_.time(1:end-1)+0.5*dtPi,dPi?./avPi?/dxPi?); epsPi?.units = ''1/km'';',1:4)

mvaE1_fromTidivNe = mvaTi1.resample(epsNe1).trace/3*units.eV/units.e*epsNe1;
c_eval('mvaE?_TiepsNe = mvaTi?.resample(epsNe?).trace/3*units.eV/units.e*epsNe?;',1:4)
c_eval('mvaE?_Pi = mvaTi?.resample(epsPi?).trace/3*units.eV/units.e*epsPi?;',1:4)
c_eval('mvaE?_Pe = mvaTe?.resample(epsPe?).trace/3*units.eV/units.e*epsPe?;',1:4)
mvaAvE_TiepsNe = (mvaE1_TiepsNe + mvaE2_TiepsNe.resample(mvaE1_TiepsNe) + mvaE3_TiepsNe.resample(mvaE1_TiepsNe) + mvaE4_TiepsNe.resample(mvaE1_TiepsNe))/4;
mvaAvE_Pi = (mvaE1_Pi + mvaE2_Pi.resample(mvaE1_Pi) + mvaE3_Pi.resample(mvaE1_Pi) + mvaE4_Pi.resample(mvaE1_Pi))/4;
mvaAvE_Pe = (mvaE1_Pe + mvaE2_Pe.resample(mvaE1_Pe) + mvaE3_Pe.resample(mvaE1_Pe) + mvaE4_Pe.resample(mvaE1_Pe))/4;

if 0 % debug mvaE?_fromPe
  %%
  figure(11)
  ic = 1;
  npanels = 7;
  h = irf_plot(npanels);
  isub = 1;
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,Pe?_);',ic)
    hca.YLabel.String = 'Pe = trace(Pe)/3';
  end
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,irf.ts_scalar(Pe?_.time(1:end-1)+0.5*dtPe,avPe?));',ic)
    hca.YLabel.String = 'av(Pe)';
  end
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,irf.ts_scalar(Pe?_.time(1:end-1)+0.5*dtPe,dPe?));',ic)
    hca.YLabel.String = 'diff(Pe)';
  end
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,irf.ts_scalar(Pe?_.time(1:end-1)+0.5*dtPe,dPe?./avPe?));',ic)
    hca.YLabel.String = 'diff(Pe)/av(Pe)';
  end
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,{irf.ts_scalar(Pe?_.time(1:end-1)+0.5*dtPe,dPe?./avPe?/dxPe?),epsPe?},''comp'');',ic)
    hca.YLabel.String = {'diff(Pe)/av(Pe)/dx','(1/km)'};
    c_eval('dt = dtPe?; dx = dxPe?;')
    irf_legend(hca,sprintf('dx = dt*v_{cs} = %g s * %g km/s = %g km',dt,CS_normal_velocity,dx),[0.05 0.98])
  end
  if 1
    hca = h(isub); isub = isub + 1;
    c_eval('irf_plot(hca,mvaE?_fromPe);',ic)
    hca.YLabel.String = 'E (mV/m)';
  end
  irf_zoom(h,'x',tintObs)
  irf_zoom(h,'y')
end
if 1
  %%
  figure(10); 
  ic = 1;
  c_eval('h=irf_plot({mvaE?.z.resample(mvaE?_fromTidivNe),mvaE?_fromPi,mvaE?_fromTidivNe,1*mvaVexB?.z,mvaOhmGradPe},''comp'');',ic)
  legend(h,sprintf('E%g',ic),sprintf('E(Pi%g,v_{CS}=%.0f km/s)',ic,CS_normal_velocity),sprintf('E(Ti%g,ne%g,v_{CS}=%.0f km/s)',ic,ic,CS_normal_velocity),sprintf('VexB%g',ic),'div(Pe)/ne');
  irf_zoom(h,'x',tintObs)
  irf_zoom(h,'y')
end
sumEOM1N = irf.ts_scalar(mvaE1_fromTidivNe.time,mvaE1.z.resample(mvaE1_fromTidivNe).data + mvaVexB1.z.resample(mvaE1_fromTidivNe).data + mvaE1_fromTidivNe.data);
%c_eval('dtP = mvaPe?.time(2)-mvaPe?.time(1); LPNN?time = mvaPe?.time(1:end-1)+0.5*dtP;',ic)
%c_eval('sumP? = mvaPe?.zz.data(1:end-1)+mvaPe?.zz.data(2:end);',ic)
%c_eval('LPNN? = irf.ts_scalar(LPNN?time,0.5*sumP?./diff(mvaPe?.zz.data(:,1))*CS_normal_velocity*dtP);',ic)

%% for paper: pressure, epsilon, En = (T/e)*epsn
ic = 3;
npanels = 9;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % av B
  hca = irf_panel('B av');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvB.x,mvaAvB.y,mvaAvB.z,mvaAvB.abs},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
end
if 1 % av n
  hca = irf_panel('n av');
  set(hca,'ColorOrder',mms_colors('12'))
  irf_plot(hca,{avNe,avNi},'comp');
  hca.YLabel.String = {'n','(cc)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % av T
  hca = irf_panel('T av');
  Tishift = 430;
  Tidiv = 10;
  set(hca,'ColorOrder',mms_colors('14b23'))
  irf_plot(hca,{avTe,facAvTe.xx,0.5*(facAvTe.yy+facAvTe.zz),avTi-Tishift,avTi/Tidiv},'comp');
  hca.YLabel.String = {'T','(ev)'};
  set(hca,'ColorOrder',mms_colors('14b23'))
  irf_legend(hca,{'T_e','T_{e||}','T_{e\perp}',sprintf('T_i - %g eV',Tishift),sprintf('T_i/%g',Tidiv)},[0.98 0.9],'fontsize',12);
end
if 1 % av P, 4sc av
  hca = irf_panel('P av');
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  irf_plot(hca,{avPe,avPi-Pishift,avPB},'comp');
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 1 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  hca.YLabel.String = {sprintf('P_%g',ic),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms1
  isc = 1;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms2
  isc = 2;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms3
  isc = 3;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % av P, 1sc, mms4
  isc = 4;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % eps ne Te Pe in z (N)
  hca = irf_panel('eps normal direction');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.z,mvaEpsTe.z,mvaEpsPe.z},'comp');
  hca.YLabel.String = {'\epsilon_N','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'\epsilon_{ne}','\epsilon_{Te}','\epsilon_{Pe}'},[0.98 0.9],'fontsize',12);
end
if 0 % eps n
  hca = irf_panel('eps n');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.x,mvaEpsNe.y,mvaEpsNe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsNi.x,mvaEpsNi.y,mvaEpsNi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla n/n','(1/km)'};
end
if 0 % eps T
  hca = irf_panel('eps T');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsTe.x,mvaEpsTe.y,mvaEpsTe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsTi.x,mvaEpsTi.y,mvaEpsTi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla T/T','(1/km)'};
end
if 0 % eps P
  hca = irf_panel('eps P');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsPe.x,mvaEpsPe.y,mvaEpsPe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsPi.x,mvaEpsPi.y,mvaEpsPi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla P/P','(1/km)'};
end
if 0 % EN, 1 sc
  hca = irf_panel('E normal direction, 2 sc');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVexB?.z,mvaOhmGradPe.z,mvaAvE_Pi.z,mvaE_Pi.z},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_exB_%.0f',ic),'div(Pe)/ne','(T_i/e)\epsilon_{ne}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN, 1 sc, ions
  hca = irf_panel('E normal direction, 1 sc, ions eom');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVixB?.z,1*mvaE?_fromPi},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_ixB_%.0f',ic),'div(Pi)/ne'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN, 1 sc, electrons
  hca = irf_panel('E normal direction, 1 sc, electrons eom');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVexB?.z,1*mvaE?_fromPe,mvaOhmGradPe.z},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_exB_%.0f',ic),'div(Pe)/ne','div(Pe)/ne (4sc)'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN av , ions
  hca = irf_panel('E ions normal direction');  
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_plot(hca,{mvaAvE.z,-1*mvaAvVixB.z,mvaEfromdivPi.z,mvaAvE_fromPi},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'4sc av','-v_ixB','(T_i/e)\epsilon_{ne}','(T_i/e)\epsilon_{Pi}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN av, electrons
  hca = irf_panel('E electrons normal direction');
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_plot(hca,{mvaE1.z,mvaE3.z,mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  %irf_plot(hca,{mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  irf_plot(hca,{mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = sprintf('MMS %g',ic);

for ip = [3 5 7]
 % h(ip).YLim = [-0.06 0.06];
end

add_length_on_top(h(1),CS_normal_velocity,0.2)

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% single spacecraft, momentum equations
ic = 3;
npanels = 9;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z,mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cc)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % av T
  hca = irf_panel('T av');
  Tishift = 430;
  Tidiv = 10;
  set(hca,'ColorOrder',mms_colors('14b23'))
  irf_plot(hca,{avTe,facAvTe.xx,0.5*(facAvTe.yy+facAvTe.zz),avTi-Tishift,avTi/Tidiv},'comp');
  hca.YLabel.String = {'T','(ev)'};
  set(hca,'ColorOrder',mms_colors('14b23'))
  irf_legend(hca,{'T_e','T_{e||}','T_{e\perp}',sprintf('T_i - %g eV',Tishift),sprintf('T_i/%g',Tidiv)},[0.98 0.9],'fontsize',12);
end
if 1 % av P, 4sc av
  hca = irf_panel('P av');
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  irf_plot(hca,{avPe,avPi-Pishift,avPB},'comp');
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 1 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  hca.YLabel.String = {sprintf('P_%g',ic),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms1
  isc = 1;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms2
  isc = 2;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % P, 1sc, mms3
  isc = 3;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % av P, 1sc, mms4
  isc = 4;
  hca = irf_panel(sprintf('P sc %g',isc));
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',isc)
  hca.YLabel.String = {sprintf('P_%g',isc),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % eps ne Te Pe in z (N)
  hca = irf_panel('eps normal direction');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.z,mvaEpsTe.z,mvaEpsPe.z},'comp');
  hca.YLabel.String = {'\epsilon_N','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'\epsilon_{ne}','\epsilon_{Te}','\epsilon_{Pe}'},[0.98 0.9],'fontsize',12);
end
if 0 % eps n
  hca = irf_panel('eps n');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.x,mvaEpsNe.y,mvaEpsNe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsNi.x,mvaEpsNi.y,mvaEpsNi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla n/n','(1/km)'};
end
if 0 % eps T
  hca = irf_panel('eps T');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsTe.x,mvaEpsTe.y,mvaEpsTe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsTi.x,mvaEpsTi.y,mvaEpsTi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla T/T','(1/km)'};
end
if 0 % eps P
  hca = irf_panel('eps P');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsPe.x,mvaEpsPe.y,mvaEpsPe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsPi.x,mvaEpsPi.y,mvaEpsPi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla P/P','(1/km)'};
end
if 0 % EN, 1 sc
  hca = irf_panel('E normal direction, 2 sc');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVexB?.z,mvaOhmGradPe.z,mvaAvE_fromPi.z},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_exB_%.0f',ic),'div(Pe)/ne','(T_i/e)\epsilon_{ne}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN, 1 sc, ions
  hca = irf_panel('E normal direction, 1 sc, ions eom');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVixB?.z,1*mvaE?_Pi,1*mvaE_TiepsPi.z,1*mvaAvE_Pi},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_ixB_%.0f',ic),'div(Pi)/ne (1 sc)','div(Pi)/ne (4sc)','<div(Pi)/ne> (4sc av)'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN, 1 sc, ions
  hca = irf_panel('E normal direction, 1sc E+vexB');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z.resample(mvaVixB?) + mvaVixB?.z,1*mvaE?_Pi,1*mvaE_TiepsPi.z,1*mvaAvE_Pi},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f+v_ixB_%.0f',ic,ic),'div(Pi)/ne (1 sc)','div(Pi)/ne (4sc)','<div(Pi)/ne> (4sc av)'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 1 % EN, 1 sc, electrons
  hca = irf_panel('E normal direction, 1 sc, electrons eom');
  %ic = 1;  
  set(hca,'ColorOrder',mms_colors('matlab'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaVexB?.z,1*mvaE?_fromPe,mvaOhmGradPe.z},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{sprintf('E_%.0f',ic),sprintf('-v_exB_%.0f',ic),'div(Pe)/ne','div(Pe)/ne (4sc)'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 0 % EN av , ions
  hca = irf_panel('E ions normal direction');  
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_plot(hca,{mvaAvE.z,-1*mvaAvVixB.z,mvaEfromdivPi.z,mvaAvE_fromPi},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'4sc av','-v_ixB','(T_i/e)\epsilon_{ne}','(T_i/e)\epsilon_{Pi}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
if 0 % EN av, electrons
  hca = irf_panel('E electrons normal direction');
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_plot(hca,{mvaE1.z,mvaE3.z,mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  %irf_plot(hca,{mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  irf_plot(hca,{mvaAvE.z,mvaEfromdivPi.z,-1*mvaAvVexB.z,mvaOhmGradPe.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms3','4sc av','(T_i/e)\epsilon_{ne}','-v_exB'},[0.98 0.9],'fontsize',12);
end
%irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = sprintf('MMS %g',ic);

for ip = [3 5 7]
 % h(ip).YLim = [-0.06 0.06];
end

add_length_on_top(h(1),CS_normal_velocity,0.2)

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Pressures
ic = 1;
npanels = 10;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 0 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % PB
  hca = irf_panel('PB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{PB1,PB2,PB3,PB4},'comp');
  hca.YLabel.String = {'P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Pi
  hca = irf_panel('Pi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaPi1.trace/3,mvaPi2.trace/3,mvaPi3.trace/3,mvaPi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Pe
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaPe1.trace/3,mvaPe2.trace/3,mvaPe3.trace/3,mvaPe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % ni
  hca = irf_panel('ni');
  set(hca,'ColorOrder',[mms_colors('1234');mms_colors('1234')])
  irf_plot(hca,{ni1.tlim(tint),ni2.tlim(tint),ni3.tlim(tint),ni4.tlim(tint),},'comp');
  hca.YLabel.String = {irf_ssub('n_{i}',ic),'(cm^{-3})'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',[mms_colors('1234');mms_colors('1234')])
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % Pressures
  hca = irf_panel('Pressure');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'mvaPi?.trace/3,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 1 % Pressures
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  if 1
  irf_plot(hca,{mvaPe1.trace/3+PB1.resample(mvaPe1)+mvaPi1.resample(mvaPe1).trace/3,...    
                mvaPe2.trace/3+PB2.resample(mvaPe2)+mvaPi1.resample(mvaPe2).trace/3,...    
                mvaPe3.trace/3+PB3.resample(mvaPe3)+mvaPi1.resample(mvaPe3).trace/3,...    
                mvaPe4.trace/3+PB4.resample(mvaPe4)+mvaPi1.resample(mvaPe4).trace/3},'comp');
  else
    irf_plot(hca,{0.5*(mvaPe1.yy+mvaPe1.yy)+PB1.resample(mvaPe1),...    
                  0.5*(mvaPe2.yy+mvaPe2.yy)+PB2.resample(mvaPe2),...    
                  0.5*(mvaPe3.yy+mvaPe3.yy)+PB3.resample(mvaPe3),...    
                  0.5*(mvaPe4.yy+mvaPe4.yy)+PB4.resample(mvaPe4)},'comp');
  end
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Pressures
  hca = irf_panel('PB + Pe');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval(['irf_plot(hca,{mvaPe?.trace/3+PB?.resample(mvaPe?)},''comp'');'],ic)
  hca.YLabel.String = {'P_e+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))  
end
if 1 % Pressures
  hca = irf_panel('PB + Pi');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval(['irf_plot(hca,{mvaPi?.trace/3+PB?.resample(mvaPi?)},''comp'');'],ic)
  hca.YLabel.String = {'P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% epsilon
ic = 1;
npanels = 7;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % av B
  hca = irf_panel('B av');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaAvB.x,mvaAvB.y,mvaAvB.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
end
if 1 % av n
  hca = irf_panel('n av');
  set(hca,'ColorOrder',mms_colors('12'))
  irf_plot(hca,{avNe,avNi},'comp');
  hca.YLabel.String = {'n','(cc)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % eps n
  hca = irf_panel('eps n');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.x,mvaEpsNe.y,mvaEpsNe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsNi.x,mvaEpsNi.y,mvaEpsNi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla n/n','(1/km)'};
end
if 1 % av T
  hca = irf_panel('T av');
  Tishift = 430;
  Tidiv = 10;
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{avTe,avTi-Tishift,avTi/Tidiv,facAvTe.xx,0.5*(facAvTe.yy+facAvTe.zz)},'comp');
  hca.YLabel.String = {'T','(ev)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'T_e',sprintf('T_i - %g eV',Tishift),sprintf('T_i/%g',Tidiv),'T_{e||}','T_{e\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % eps T
  hca = irf_panel('eps T');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsTe.x,mvaEpsTe.y,mvaEpsTe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsTi.x,mvaEpsTi.y,mvaEpsTi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla T/T','(1/km)'};
end
if 1 % av P
  hca = irf_panel('P av');
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  irf_plot(hca,{avPe,avPi-Pishift,avPB},'comp');
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 1 % eps P
  hca = irf_panel('eps P');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsPe.x,mvaEpsPe.y,mvaEpsPe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsPi.x,mvaEpsPi.y,mvaEpsPi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla P/P','(1/km)'};
end
irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

for ip = [3 5 7]
 % h(ip).YLim = [-0.06 0.06];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% energy selection, electrostatic potential
ic = 1;
npanels = 5;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

tref = irf_time('2015-11-12T07:19:21.60Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)
c_eval('mvaPhi?N = mvaIntE?N*CS_normal_velocity*-1;',ic)


scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % av B
  hca = irf_panel('B av');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaAvB.x,mvaAvB.y,mvaAvB.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
end
if 0 % av n
  hca = irf_panel('n av');
  set(hca,'ColorOrder',mms_colors('12'))
  irf_plot(hca,{avNe,avNi},'comp');
  hca.YLabel.String = {'n','(cc)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 0 % eps n
  hca = irf_panel('eps n');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsNe.x,mvaEpsNe.y,mvaEpsNe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsNi.x,mvaEpsNi.y,mvaEpsNi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla n/n','(1/km)'};
end
if 1 % av T
  hca = irf_panel('T av');
  Tishift = 430;
  Tidiv = 10;
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{avTe,avTi-Tishift,avTi/Tidiv,facAvTe.xx,0.5*(facAvTe.yy+facAvTe.zz)},'comp');
  hca.YLabel.String = {'T','(ev)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'T_e',sprintf('T_i - %g eV',Tishift),sprintf('T_i/%g',Tidiv),'T_{e||}','T_{e\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % eps T
  hca = irf_panel('eps T');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsTe.x,mvaEpsTe.y,mvaEpsTe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsTi.x,mvaEpsTi.y,mvaEpsTi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla T/T','(1/km)'};
end
if 0 % av P
  hca = irf_panel('P av');
  Pishift = 0.42;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  irf_plot(hca,{avPe,avPi-Pishift,avPB},'comp');
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
end
if 0 % eps P
  hca = irf_panel('eps P');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaEpsPe.x,mvaEpsPe.y,mvaEpsPe.z},'comp');
  hold(hca,'on')
  irf_plot(hca,{mvaEpsPi.x,mvaEpsPi.y,mvaEpsPi.z},'comp','--');
  hold(hca,'off')
  hca.YLabel.String = {'\nabla P/P','(1/km)'};
end
if 0 % E
  hca = irf_panel('E av');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaAvE.x.tlim(tint),mvaAvE.y.tlim(tint),mvaAvE.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 1 % potential in normal direction
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('z'))
  c_eval('irf_plot(hca,{mvaPhi?N},''comp'');',ic)
  hca.YLabel.String = {'\phi','(V)'};
  set(hca,'ColorOrder',mms_colors('z'))
  irf_legend(hca,{'N'},[0.98 0.9],'fontsize',12);  
  irf_pl_mark(hca,tref,'k')
end
irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

for ip = [3 5 7]
 % h(ip).YLim = [-0.06 0.06];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end