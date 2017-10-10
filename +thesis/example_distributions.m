%% Load data





%% Plot distributions


%time1 = irf_time('2015-12-06T23:38:31.588Z','utc>epochtt'); % Yuri's field aligned beam
%time2 = irf_time('2015-10-16T10:33:45.238Z','utc>epochtt'); % Dougnut close to fluxrope
%time3 = irf_time('2015-10-16T10:33:30.326Z','utc>epochtt'); % Outflow crescent
%time4 = irf_time('2015-10-16T10:33:30.326Z','utc>epochtt'); % Flat top
%time5 = irf_time('2015-10-16T10:33:30.326Z','utc>epochtt'); % Loss cone
events = {};
if 1 % Adiabatic electron inflow in Nov 12 event.
  event.sc = 4;
  event.time = irf_time('2015-11-12T07:19:21.342Z','utc>epochtt');
  event.tint = event.time + [-10 10];
  events{end+1} = event;
end
if 0 % Bulk drift in Nov 12 event.
  event.sc = 4;
  event.time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');
  event.tint = event.time + [-10 10];
  events{end+1} = event;
end
if 0 % flat top
  event.sc = 1;
  event.time = irf_time('2015-11-12T07:19:17.70Z','utc>epochtt');
  event.tint = event.time + [-10 10];
  events{end+1} = event;
end
if 0 % Yuri's field aligned beam
  event.sc = 4;
  event.time = irf_time('2015-12-06T23:38:31.588Z','utc>epochtt');
  event.tint = event.time + [-10 10];
  events{end+1} = event;
end

nRows = 1;
nCols = 1;
for isub = 1:nRows*nCols; h(isub) =  subplot(nRows,nCols,isub); end

ic = 4;
loadData = 0;
calculatePitchangles = 0;
isub = 1;
for iTime = 1:numel(events);
  event = events{iTime};
  time = event.time;
  tint = event.tint;
  ic = event.sc;
  if 1
    hca = h(isub); isub = isub + 1;    
    if loadData
      c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
      c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
      c_eval('dmpaB?slow = dmpaB?.resample(ePDist?);',ic)
    end
    if calculatePitchangles    
      c_eval('ePitch?par = ePDist?.pitchangles(dmpaB?,[0 15]).convertto(''s^3/km^6'');',ic)
      c_eval('ePitch?perp = ePDist?.pitchangles(dmpaB?,[75 105]).convertto(''s^3/km^6'');',ic)
      c_eval('ePitch?apar = ePDist?.pitchangles(dmpaB?,[165 180]).convertto(''s^3/km^6'');',ic)
    end
   
    n = [7];
    m = [0];
    t = [60];
    vd = [0];
    d = [1];
    a1 = [1];
    a2 = [0];
    toPlot = [1];
    
    %f = cn.maxwellian(ve,35,20,0,'e','3D'); % f = cn.maxwellian(v,T,n,vd,species,optional);      
    
    %v = linspace(cn_eV2v(10,'eV'),cn_eV2v(5e3,'eV'),100);
    %plot(hca,cn.maxwellian(v,600,20,0,'e','3D')*1e18);
    
    c_eval('it = find(abs(ePDist?.time-time)==min(abs(ePDist?.time-time)));',ic);
    c_eval('time = ePDist?(it).time;',ic)
    timeUTC = time.utc;
    c_eval('hl = plot(hca,ePitch?par(it).depend{1},[ePitch?par(it).data; ePitch?perp(it).data; ePitch?apar(it).data])',ic)
    
    hca.YScale = 'log'; hca.XScale = 'log';
    
    hca.YLim = [1e0 1e5]; hca.XLim = [1e1 1e3];
    hca.YTick = [1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
    %legend(hca,num2str(ePitch3.depend{2}(indPA)'),'location','southwest')
    
    hold(hca,'on')
    hh=whamp.plot_f(hca,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
    hh.Children(1).Color = [0 0 0];
    hold(hca,'off')
    hca.YGrid = 'off';
    hca.XGrid = 'off';
    hca.Title.String = [timeUTC(1:10) '  ' timeUTC(12:23)];
    hca.YLabel.String = 'f_e (s^3/km^6)'; hca.XLabel.String = 'E (eV)';
  end
end
for  ii = 4:6; hh.Children(ii).LineWidth = 2; end
hh.Children(1).LineStyle = '--'; 
hh.Children(2).Visible = 'off'; 
hh.Children(3).Visible = 'off';
irf_legend(hca,'-- Maxwellian',[0.1 0.2],'fontsize',16)
irf_legend(hca,sprintf('T_e = %.0f eV, n_e = %.0f cm^{-3}',t,n),[0.1 0.1],'fontsize',16)
for isub = 1:nRows*nCols
  hca = h(isub);
  hca.FontSize = 14;
end
hca.Position