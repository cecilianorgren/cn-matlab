%% Le20?? equations of state
tref = irf_time('2017-07-11T22:33:25.00Z','utc>epochtt');
c_eval('cglPe?par = pi/6*(ne?/ne?.resample(tref)).^3.*(gseB?.abs.resample(tref).data/gseB?.resample(ne?).abs).^2;',ic)
%c_eval('Pr?par'
%cglPe1par = 
%aa = pi/6*(ne1/ne1.resample(tref)).^3;
%bb = (gseB1.abs.resample(tref).data/gseB1.resample(ne1).abs).^2;

%% Le2010, Magnitude of the Hall fields during magnetic reconnection
tref1 = irf_time('2017-07-11T22:33:21.00Z','utc>epochtt'); % left of cs
tref2 = irf_time('2017-07-11T22:33:27.00Z','utc>epochtt'); % right of cs
tref3 = irf_time('2017-07-11T22:33:32.00Z','utc>epochtt'); % left of cs
iref = 1:3;

%c_eval('eBetainf! = 2*gsePe?.resample(tref!).trace.data/3/gseB?.abs2.resample(tref!).data;',ic,iref)
c_eval('eBetainf! = beta?e.resample(tref!).data;',ic,iref)
c_eval('Binf! = gseB?.abs.resample(tref!).data;',ic,iref)
c_eval('nhat! = ne?/ne?.resample(tref!).data; nhat!.name = ''nhat (Le2010)'';',ic,iref)
Bhref = 12;
c_eval('BH! = Bhref-Binf!*(pi*nhat!.^3*eBetainf!/12).^(1/4); BH!.name = ''B_H (Le2010)''; ?;',ic,iref)

tintZoom = irf.tint('2017-07-11T22:33:15.00Z/2017-07-11T22:33:45.00Z');


ic = 1;
npanels = 3;
h = irf_plot(npanels);

it = 1;

zoomy = [];
isub = 0;
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
end
if 0 % nhat, one, defined by 'it'
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('nhat');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{nhat?},''comp'')',ic); 
  ylabel(hca,'n/n_{ref}','interpreter','tex');  
end
if 1 % nhat
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('nhat');
  set(hca,'ColorOrder',mms_colors('123'))
  irf_plot(hca,{nhat1,nhat2,nhat3},'comp');
  ylabel(hca,'n/n_{ref}','interpreter','tex'); 
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'n_{ref,1}','n_{ref,2}','n_{ref,3}'},[0.98 0.99],'fontsize',12);  
end
if 0 % B, one, defined by 'it'  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B M comp');
  set(hca,'ColorOrder',mms_colors('y1'))
  c_eval('irf_plot(hca,{mvaB?.y,BH1},''comp'');',ic)
  hca.YLabel.String = {'B_H','(nT)'};
  set(hca,'ColorOrder',mms_colors('y1'))
  irf_legend(hca,{'B_M^{obs}','B_H (Le2010)'},[0.98 0.9],'fontsize',12);  
end
if 1 % B 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B M comp');
  set(hca,'ColorOrder',mms_colors('123b'))
  c_eval('irf_plot(hca,{BH1,BH2,BH3,mvaB?.y},''comp'');',ic)
  hca.YLabel.String = {'B_H','(nT)'};
  set(hca,'ColorOrder',mms_colors('123b'))
  irf_legend(hca,{'B_{H,ref1}','B_{H,ref2}','B_{H,ref3}','B_M^{obs}'},[0.98 0.99],'fontsize',12);  
  irf_legend(hca,{sprintf('%g nT-B_{H,Le2010}',Bhref)},[0.01 0.9],'fontsize',12,'color',[0 0 0]);  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
c_eval('hmark = irf_pl_mark(h(?),tref!.epochUnix,mms_colors(''!''));',1:npanels,iref)


%% Le2009/Le2010, CGL like scalings, BH scaling
tref1 = irf_time('2017-07-11T22:33:27.00Z','utc>epochtt'); % left of cs
tref2 = irf_time('2017-07-11T22:33:31.00Z','utc>epochtt'); % right of cs
c_eval('nhat?_! = ne?/ne?.resample(tref!).data; nhat?.name = ''nhat (Le2010)''; nhat?.userData.LABLAXIS = ''nhat'';',ic,1:2)
c_eval('Bhat?_! = gseB?.abs/gseB?.abs.resample(tref!).data; Bhat?.name = ''Bhat'';',ic,1:2)
c_eval('pref?_! = facPe?.trace.resample(tref!).data/3;',ic,1:2)
c_eval('cglPe?par_! =  pref?_!*pi/6*nhat?_!.^3/(Bhat?_!.resample(nhat?_!).abs.^2); cglPe?par.name = ''cgl Pe par (Le2010)'';',ic,1:2)
% version below has no pi/6 term, which makes it match up better with data
c_eval('cglPe?par_! =  pref?_!*nhat?_!.^3/(Bhat?_!.resample(nhat?_!).abs.^2); cglPe?par.name = ''cgl Pe par (Le2010)'';',ic,1:2)
c_eval('cglPe?perp_! = pref?_!*nhat?_!*Bhat?_!.resample(nhat?_!).abs; cglPe?perp.name = ''cgl Pe perp (Le2010)'';',ic,1:2)
c_eval('cglPe?parperp_! = cglPe?par_!/cglPe?perp_!; cglPe?parperp.name = ''cgl Pe par/perp (Le2010)'';',ic,1:2)

c_eval('Pdiff? = irf.ts_scalar(mvaPe?.time,units.mu0*1e9*(mvaPe?.xx.data-0.5*(mvaPe?.yy.data+mvaPe?.zz.data))); Pdiff?.name = ''Ppar-Pperp''; Pdiff?.units = ''(nT)^2'';',ic)
c_eval('tens? = gseB?.abs2.resample(Pdiff?)-Pdiff?;',ic)


ic = 1;
npanels = 9;
h = irf_plot(npanels);

zoomy = [];
isub = 0;
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Bhat
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Bhat');
	set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{Bhat?_1,Bhat?_2},''comp'')',ic); 
  ylabel(hca,'B/B_{ref}','interpreter','tex');
	set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'B_{ref,1}','B_{ref,2}'},[0.98 0.9],'fontsize',12);  
end
if 1 % nhat
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('nhat');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{nhat?_1,nhat?_2},''comp'')',ic); 
  ylabel(hca,'n/n_{ref}','interpreter','tex');
  irf_legend(hca,{'n_{ref,1}','n_{ref,2}'},[0.98 0.9],'fontsize',12);  
end
if 1
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('P perp');
  set(hca,'ColorOrder',mms_colors('x12'))
  c_eval('irf_plot(hca,{0.5*(mvaPe?.yy+mvaPe?.zz),cglPe?perp_1,cglPe?perp_2},''comp'');',ic); 
  ylabel(hca,{'P_{\perp}','(nPa)'},'interpreter','tex'); 
  set(hca,'ColorOrder',mms_colors('x12'))
  irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])
end
if 1
  isub = isub + 1;
  zoomy = [zoomy isub];    
  hca = irf_panel('P par');
  set(hca,'ColorOrder',mms_colors('x12'))
  c_eval('irf_plot(hca,{mvaPe?.xx,cglPe?par_1,cglPe?par_2},''comp'');',ic); 
  ylabel(hca,{'P_{||}','(nPa)'},'interpreter','tex'); 
  set(hca,'ColorOrder',mms_colors('x12'))
  irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])
end
if 1
  isub = isub + 1;
  zoomy = [zoomy isub];    
  set(hca,'ColorOrder',mms_colors('x12'))
  hca = irf_panel('P par/perp');
  c_eval('irf_plot(hca,{mvaPe?.xx/(0.5*(mvaPe?.yy+mvaPe?.zz)),cglPe?parperp_1,cglPe?parperp_2},''comp'');',ic); 
  ylabel(hca,{'P_{||}/P_{\perp}',''},'interpreter','tex'); 
  set(hca,'ColorOrder',mms_colors('x12'))
  irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])
end
if 1
  isub = isub + 1;
  zoomy = [zoomy isub];   
  set(hca,'ColorOrder',mms_colors('3'))
  hca = irf_panel('mu0*(Ppar+PperpP)');
  c_eval('irf_plot(hca,Pdiff?);',ic); 
  ylabel(hca,{'\mu_0(P_{||}-P_{\perp})','(nT)'},'interpreter','tex'); 
end
if 1
  isub = isub + 1;
  zoomy = [zoomy isub];   
  set(hca,'ColorOrder',mms_colors('3'))
  hca = irf_panel('B-Ppar+PperpP');
  c_eval('irf_plot(hca,tens?);',ic); 
  ylabel(hca,{'B^2-\mu_0(P_{||}-P_{\perp})','(nT)'},'interpreter','tex'); 
end
if 0
  isub = isub + 1;
  zoomy = [zoomy isub];   
  hca = irf_panel('curv B');
  c_eval('irf_plot(hca,mvaCurvB);',ic); 
  ylabel(hca,{'curv b)','(1/km)'},'interpreter','tex'); 
end

c_eval('hmark1 = irf_pl_mark(h(?),tref1.epochUnix,mms_colors(''1''));',1:npanels)
c_eval('hmark2 = irf_pl_mark(h(?),tref2.epochUnix,mms_colors(''2''));',1:npanels)

irf_zoom(h,'x',colon(tref1,tref2-tref1,tref2)+0.5*[-1 1])
irf_zoom(h(zoomy),'y')
%hca = irf_panel('P par'); hca.YLim = [0 0.2];
%hca = irf_panel('P par/perp'); hca.YLim = [0 4];