%tDRCenter = irf.tint('2015-10-16T10:33:28.00Z',0.00001); tDRCenter = tDRCenter(1);
%tint = tDRCenter+0.06*[-1 1] ;
%tint = irf.tint('2015-10-16T10:33:30.23Z',0.04);
%mms.plot_skymap(gca,desDist1,'tint',tint,'energy',80);

timeOverview = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:35.00Z');
timeIntervals = {irf.tint('2015-10-16T10:33:30.15Z',0.033),...
                 irf.tint('2015-10-16T10:33:30.25Z',0.033),...
                 irf.tint('2015-10-16T10:33:30.35Z',0.033),...
                 irf.tint('2015-10-16T10:33:32.00Z',0.033),...
                };
              
timeIntervals = {irf.tint('2015-10-16T10:33:26.00Z',0.033),...
                 irf.tint('2015-10-16T10:33:26.50Z',0.033),...
                 irf.tint('2015-10-16T10:33:26.90Z',0.033),...
                 irf.tint('2015-10-16T10:33:27.50Z',0.033),...
                };              
%%
toplot = 5;

%% Set up plot
%[h1,h2] = mms.initialize_comined_plot(8,4,4,5,'horizontal');
[h1,h2] = mms.initialize_comined_plot(8,4,4,5,'vertical');

%% Plot time series
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');

ic = 1;

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('brst E');
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E_{DSL}','(mV/m)'};
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

if 0
hca = irf_panel('brst scPot');
c_eval('irf_plot(hca,(-1)*P?brst.tlim(tint));',ic);
hca.YLabel.String = {'-scPot','(V)'};
end

hca = irf_panel('brst n');
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
c_eval('irf_plot(hca,ne?_lowres.tlim(tint));',ic);
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
irf_legend(hca,{'n_e','n_i'},[0.95 0.95]);

hca = irf_panel('brst Te');
c_eval('irf_plot(hca,Te?_lowres.tlim(tint));',ic);
hca.YLabel.String = {'T','(eV)'};
hca.YScale = 'lin';
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.95 0.95]);

hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst.tlim(tint).tlim(tint));',ic);
hca.YLabel.String = {'v_e','(km/s)'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
hold(hca,'on');
c_eval('irf_plot(hca,Te?_lowres.tlim(tint))',ic)
%c_eval('irf_plot(hca,veeV?)',ic)
hold(hca,'off');
irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.5 0.1])
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%caxis(hca,[-15 -5])
hca.YLim = [10 30e3];
ylabel(hca,{'E_{e,OMNI}','(eV)'},'fontsize',12,'Interpreter','tex');

% Electron velocity in eV
me = 9.1094e-31;
e = 1.6022e-19;
c_eval('veeV? = ve?brst.abs*ve?brst.abs*1e6*0.5*me/e;',ic);
c_eval('Bhat = dmpaB?brst.resample(veeV?.tlim(tint).time); Bhat = Bhat/Bhat.abs',ic);
c_eval('veeV?par = ve?brst.dot(Bhat).times(ve?brst.dot(Bhat))*1e6*0.5*me/e;',ic)
hca=irf_panel('edistparperp');
c_eval('irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
hold(hca,'on')
c_eval('irf_plot(hca,veeV?par)',ic)
hold(hca,'off')
irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,1*[-1 1])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E','(eV)'});
colormap(hca,cn.cmap('bluered3'))

hca=irf_panel('edistparapar');
c_eval('irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,1*[-1 1])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E','(eV)'});
colormap(hca,cn.cmap('bluered3'))


h1(1).Title.String = ['MMS ' num2str(ic)];
irf_zoom(h1,'x',tint)
irf_plot_axis_align
irf_zoom(h1(1:5),'y')
for ii = 6:8
  h1(ii).YLim = [10 500];
end

irf_zoom(h1,'x',timeOverview)
irf_plot_axis_align
irf_zoom(h1(1:5),'y')  
hmark = 1:5;
%% Plot particle skymaps
% 1: skymaps and psds for 3 time intervals, 'first En'
% 2: skymaps and psds for 3 time intervals, 'firstfirst En'
% 3: psd for 6 time intervals, 'first En'
% 4:
% 5:
% 6: 
% 7:
% 8:

if 0 % get mva direction
  %%
  % Automatic results made from gui above
  ic = 1:4;
  tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
  c_eval('[out?,l?,mva?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
  mva = mva3;
  % rotate around x-axis 
  turn_angle = 0; % degrees
  tmpv = mva;
  mva(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
  mva(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
  % v is 3x3 matrix:
  %v1: [0.0906 0.0896 0.9918] - first row
  %v2: [0.3510 -0.9349 0.0524] - second row
  %v3: [0.9320 0.3434 -0.1161] - third row   
end

switch toplot
  case 1
    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:35.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y')  
    
    hmark = 1:5;
    isub = 1;
    skymapEnergy = 230;
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.10Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));    
    vectors = {hatB0,'B'};
    if norm(E0)>2, vectors(end+1,1:2) = {hatE0,'E'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.30Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));
    vectors = {hatB0,'B'};
    if norm(E0)>2, vectors(end+1,1:2) = {hatE0,'E'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.40Z',0.04);   
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));
    vectors = {hatB0,'B'};
    if norm(E0)>2, vectors(end+1,1:2) = {hatE0,'E'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
          
  case 2
    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:35.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y')  
    
    hmark = 1:5;
    isub = 1;
    skymapEnergy = 200;
    limE = 2;
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:26.10Z',0.04);   
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0)); hatE0(3) = 0;
    vectors = {hatB0,'B'};
    if norm(E0)>limE, vectors(end+1,1:2) = {hatE0,'E_{z=0}'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:26.40Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));    
    vectors = {hatB0,'B'};
    if norm(E0)>limE, vectors(end+1,1:2) = {hatE0,'E'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:27.28Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',skymapEnergy);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));
    vectors = {hatB0,'B'};
    if norm(E0)>limE, vectors(end+1,1:2) = {hatE0,'E_{z=0}'}; end
    c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
  case 3

    isub = 1;  

    if 0
      tint = irf.tint('2015-10-16T10:33:24.00Z',0.04);
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 0
      tint = irf.tint('2015-10-16T10:33:26.70Z',0.04);
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 0
      tint = irf.tint('2015-10-16T10:33:30.00Z',0.04);
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 0
      tint = irf.tint('2015-10-16T10:33:28.70Z',0.04);
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
  
    tint = irf.tint('2015-10-16T10:33:30.10Z',0.04);
    hca = h2(isub); isub = isub + 1;
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')

    tint = irf.tint('2015-10-16T10:33:30.20Z',0.04);
    hca = h2(isub); isub = isub + 1;
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')

    tint = irf.tint('2015-10-16T10:33:30.30Z',0.04);
    hca = h2(isub); isub = isub + 1; 
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P1);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')

    tint = irf.tint('2015-10-16T10:33:30.40Z',0.04);
    hca = h2(isub); isub = isub + 1; 
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P1);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')

    if 1
    tint = irf.tint('2015-10-16T10:33:30.50Z',0.04);
    hca = h2(isub); isub = isub + 1; 
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P1);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    
    tint = irf.tint('2015-10-16T10:33:30.70Z',0.04);
    hca = h2(isub); isub = isub + 1; 
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P1);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    for ii = 1:6
      h2(ii).YLim = [1e-3 2e5];
      h2(ii).YTick = 10.^[-4:2:6];
      h2(ii).XTick = 10.^[1:1:5];
    end 

    %irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:35.00Z'))
    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:27.00Z/2015-10-16T10:33:33.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y') 
  
  case 4
    isub = 1;  

    hmark = 1:5;

    if 1
      tint = irf.tint('2015-10-16T10:33:44.50Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 1
      tint = irf.tint('2015-10-16T10:33:44.75Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 1
      tint = irf.tint('2015-10-16T10:33:45.00Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 1
      tint = irf.tint('2015-10-16T10:33:45.25Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 1
      tint = irf.tint('2015-10-16T10:33:45.50Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end
    if 1
      tint = irf.tint('2015-10-16T10:33:45.75Z',0.04);
      hca = h2(isub); isub = isub + 1; 
      mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
    end

    for ii = 1:6
      h2(ii).YLim = [1e-3 2e5];
      h2(ii).YTick = 10.^[-4:2:6];
      h2(ii).XTick = 10.^[1:1:5];
    end 

    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:35.00Z/2015-10-16T10:33:50.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y') 
 
  case 5
    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:35.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y')  
    
    hmark = 1:5;
    isub = 1;
    skymapEnergy = 200;
    limE = 10;
    vlim =10*1e3; % km/s
    
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.15Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0)); hatE0(3) = 0;    
    if norm(E0)>limE, vectors = {hatE0,'E_{z=0}'}; else vectors = {}; end
    c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',[-1 0 0; 0 -1 0;hatB0],''clim'',[0 5],''vlim'',vlim,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.25Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));    
    vectors = {hatB0,'B'};
    if norm(E0)>limE, vectors = {hatE0,'E_{z=0}'}; else vectors = {}; end
    c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',[-1 0 0; 0 -1 0;hatB0],''clim'',[0 5],''vlim'',vlim,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
    
    hca = h2(isub); isub = isub + 1;
    tint = irf.tint('2015-10-16T10:33:30.35Z',0.04);    
    c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?);',ic)
    hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    hca = h2(isub); isub = isub + 1;
    c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0));
    vectors = {hatB0,'B'};
    if norm(E0)>limE, vectors = {hatE0,'E_{z=0}'}; else vectors = {[1 0 0],'X';[0 1 0],'Y';[0 0 1],'Z'}; end
    c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',[-1 0 0; 0 -1 0;hatB0],''clim'',[0 5],''vlim'',vlim,''vectors'',vectors);',ic)
    %c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',[1 0 0; 0 1 0;0 0 1],''clim'',[0 5],''vlim'',vlim,''vectors'',vectors);',ic)
    irf_pl_mark(h1(hmark),tint.epochUnix','green')
  
  case 6
    irf_zoom(h1,'x',irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:35.00Z'))
    irf_plot_axis_align
    irf_zoom(h1(1:5),'y')  
    
    hmark = 1:5;
    isub = 1;
    cs_ylim = [1e-3 2e5];
    skymapEnergy = 90;
    limE = 3;
    vlim =10*1e3; % km/s
    lmn = mva;
    mnl = [mva(2,:); mva(3,:); mva(1,:)];
    nlm = [mva(3 ,:); mva(1,:); mva(2,:)];  
    ellim = 20;
    
    for it = 1:numel(timeIntervals)  
      vectorsXYZ = {[1 0 0],'X';[0 1 0],'Y';[0 0 1],'Z'};
      vectorsLMN = {mva(1,:),'L';mva(2,:),'M';mva(3,:),'N'};
      tint = timeIntervals{it};
      c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); hatB0 = double(irf_norm(B0));
      c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); hatE0 = double(irf_norm(E0)); hatE0(3) = 0;   
      vectorsXYZ(end+1,:) = {hatB0,'B'}; vectorsLMN(end+1,:) = {hatB0,'B'};
      if norm(E0)>limE, vectorsXYZ(end+1,1:2) = {hatE0,'E'}; vectorsLMN(end+1,1:2) = {hatE0,'E'}; end      
      
      
      % Plot psd 0 90 180?
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_cross_section_psd(hca,desDist?,dmpaB?brst,''tint'',tint,''scPot'',P?,''ylim'',cs_ylim);',ic)
      hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      %hca.Title.String = [];
      
      % Plot projection onto plane
      hca = h2(isub); isub = isub + 1;
      % LM-plane
      %if norm(E0)>limE, vectors = {hatE0,'E_{z=0}'}; else vectors = {}; end
      c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',lmn,''clim'',[0 5],''vlim'',vlim,''vectors'',vectorsXYZ,''vlabel'',{''v_L'',''v_M'',''v_N''},''elevationlim'',ellim);',ic)
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
      hca.Title.String = [];
      % 
      hca = h2(isub); isub = isub + 1;
      % NL-plane
      %if norm(E0)>limE, vectors = {hatE0,'E_{z=0}'}; else vectors = {}; end
      c_eval('mms.plot_projection(hca,desDist?,''tint'',tint,''xyz'',nlm,''clim'',[0 5],''vlim'',vlim,''vectors'',vectorsXYZ,''vlabel'',{''v_N'',''v_L'',''v_M''},''elevationlim'',ellim);',ic)
      irf_pl_mark(h1(hmark),tint.epochUnix','green')
      hca.Title.String = [];
      
      % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;      
      c_eval('mms.plot_skymap(hca,desDist?,''tint'',tint,''energy'',skymapEnergy,''vectors'',vectorsLMN);',ic)
      hca.Title.String = hca.Title.String{2};
      
      % Mark time interval in time series plot
      irf_pl_mark(h1(hmark),tint.epochUnix','green')  
    end
    
    
end