%% Le20?? equations of state
tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
tref = irf_time('2015-11-12T07:19:30.00Z','utc>epochtt');
c_eval('cglPe?par = pi/6*(ne?/ne?.resample(tref)).^3.*(gseB?.abs.resample(tref).data/gseB?.resample(ne?).abs).^2;',ic)
%c_eval('Pr?par'
%cglPe1par = 
%aa = pi/6*(ne1/ne1.resample(tref)).^3;
%bb = (gseB1.abs.resample(tref).data/gseB1.resample(ne1).abs).^2;

%% Le2010, Magnitude of the Hall fields during magnetic reconnection
tref = irf_time('2015-11-12T07:19:20.00Z','utc>epochtt'); % left of cs
%tref = irf_time('2015-11-12T07:19:22.50Z','utc>epochtt'); % right of cs
c_eval('eBetainf? = 2*gsePe?.resample(tref).trace.data/3/gseB?.abs.resample(tref).data;',ic)
c_eval('Binf? = gseB?.abs.resample(tref).data;',ic)
c_eval('nhat? = ne?/ne?.resample(tref).data; nhat?.name = ''nhat (Le2010)'';',ic)
c_eval('BH? = 2*Binf?*(pi*nhat?.^3*eBetainf?/12).^(1/4); BH?.name = ''B_H (Le2010)'';',ic)

irf_plot({mvaB1,mvaB1.y,BH1,nhat1})

%% Le2009/Le2010, CGL like scalings, BH scaling
tref1 = irf_time('2015-11-12T07:19:20.00Z','utc>epochtt'); % left of cs
tref2 = irf_time('2015-11-12T07:19:22.50Z','utc>epochtt'); % right of cs
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
h = irf_plot(10);
hca = irf_panel('B');
c_eval('irf_plot(hca,mvaB?)',ic); 
ylabel(hca,{'B','(nT)'},'interpreter','tex');

hca = irf_panel('n');
c_eval('irf_plot(hca,ne?)',ic); 
ylabel(hca,{'n','(cc)'},'interpreter','tex');

hca = irf_panel('Bhat');
c_eval('irf_plot(hca,{Bhat?_1,Bhat?_2},''comp'')',ic); 
ylabel(hca,'B/B_{ref}','interpreter','tex');

hca = irf_panel('nhat');
c_eval('irf_plot(hca,{nhat?_1,nhat?_2},''comp'')',ic); 
ylabel(hca,'n/n_{ref}','interpreter','tex');

hca = irf_panel('P perp');
c_eval('irf_plot(hca,{0.5*(mvaPe?.yy+mvaPe?.zz),cglPe?perp_1,cglPe?perp_2},''comp'');',ic); 
ylabel(hca,{'P_{\perp}','(nPa)'},'interpreter','tex'); 
irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])

hca = irf_panel('P par');
c_eval('irf_plot(hca,{mvaPe?.xx,cglPe?par_1,cglPe?par_2},''comp'');',ic); 
ylabel(hca,{'P_{||}','(nPa)'},'interpreter','tex'); 
irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])

hca = irf_panel('P par/perp');
c_eval('irf_plot(hca,{mvaPe?.xx/(0.5*(mvaPe?.yy+mvaPe?.zz)),cglPe?parperp_1,cglPe?parperp_2},''comp'');',ic); 
ylabel(hca,{'P_{||}/P_{\perp}',''},'interpreter','tex'); 
irf_legend(hca,{'FPI','CGL ref1','CGL ref2'},[0.01,0.98])

hca = irf_panel('mu0*(Ppar+PperpP)');
c_eval('irf_plot(hca,Pdiff?);',ic); 
ylabel(hca,{'\mu_0(P_{||}-P_{\perp})','(nT)'},'interpreter','tex'); 

hca = irf_panel('B-Ppar+PperpP');
c_eval('irf_plot(hca,tens?);',ic); 
ylabel(hca,{'B^2-\mu_0(P_{||}-P_{\perp})','(nT)'},'interpreter','tex'); 

hca = irf_panel('curv B');
c_eval('irf_plot(hca,mvaCurvB);',ic); 
ylabel(hca,{'curv b)','(1/km)'},'interpreter','tex'); 


hmark1 = irf_pl_mark(h,tref1.epochUnix,hca.ColorOrder(2,:));
hmark2 = irf_pl_mark(h,tref2.epochUnix,hca.ColorOrder(3,:));

irf_zoom(h,'x',colon(tref1,tref2-tref1,tref2)+0.5*[-1 1])
irf_zoom(h,'y')
hca = irf_panel('P par'); hca.YLim = [0 0.2];
hca = irf_panel('P par/perp'); hca.YLim = [0 4];