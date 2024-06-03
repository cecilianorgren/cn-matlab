%% Set up database
% units = irf_units;
irf.log('critical')
ic = 2;

%localuser = 'cno062';
localuser = 'cecilia';
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
%mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
mms.db_init('local_file_db',['/Volumes/mms']);
db_info = datastore('mms_db');   

%% Specify events and load data
% Specify events
% save('/Users/cecilia/Data/MMS/Matlab/bifurcation.mat')
clear events
iEvent = 0;

if 1 % Oieroset 2021  
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Oierøset et al. (2021)';
  events(iEvent).tint = irf.tint('2017-10-05T03:55:16.00Z/2017-10-05T03:55:20.00Z');
  events(iEvent).t_center = irf_time('2017-10-05T03:55:18.00Z','utc>EpochTT');
  events(iEvent).dt_center = [0.6000    0.0500   -0.3500   -0.8200];
  events(iEvent).v_norm = 220;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [-0.0003, 0.997, 0.0632; 0.774, -0.0395, 0.631; 0.632, 0.0486, -0.772];  
end
if 1 % Norgren 2018
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Norgren et al. (2018)';
  events(iEvent).tint = irf.tint('2015-11-12T07:19:20.50Z/2015-11-12T07:19:22.00Z');
  events(iEvent).t_center = irf_time('2015-11-12T07:19:21.20Z','utc>EpochTT');
  events(iEvent).dt_center = [0 0.05 -0.1 -0.08];
  events(iEvent).v_norm = 70;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [0.14, -0.71, 0.69; -0.54, -0.62, -0.54; 0.81, -0.29, -0.47];  
end
if 1 % Li 2019 
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Li et al. (2019)';
  events(iEvent).tint = irf.tint('2017-08-10T12:18:20.00Z/2017-08-10T12:18:45.00Z');
  events(iEvent).t_center = irf_time('2017-10-05T03:55:18.00Z','utc>EpochTT');
  events(iEvent).dt_center = [0.6000    0.0500   -0.3500   -0.8200];
  events(iEvent).v_norm = 220;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [-0.0003, 0.997, 0.0632; 0.774, -0.0395, 0.631; 0.632, 0.0486, -0.772];  
end

if 1 % Wang 2020 - Event 1
  iEvent = iEvent + 1;
  events(iEvent).paper = {'Wang et al. (2020)','     Event I        '};
  events(iEvent).tint = irf.tint('2017-06-17T20:24:03.00Z/2017-06-17T20:24:11.00Z');
  events(iEvent).t_center = irf_time('2017-06-17T20:24:07.00Z','utc>EpochTT');
  events(iEvent).dt_center = [0.4 0.08 0 0.08];
  events(iEvent).v_norm = 67;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [0.9477, 0.3023, -0.1029; -0.0855, -0.0703, -0.9939; -0.3076, 0.9506, -0.0408];    
end
if 1 % Zhou 2019 (also Wang 2020 - Event 2, who has slightly different numbers)
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Zhou et al. (2019)';
  events(iEvent).tint = irf.tint('2017-08-10T12:18:25.00Z/2017-08-10T12:18:40.00Z');
  events(iEvent).t_center = irf_time('2017-08-10T12:18:33.00Z','utc>EpochTT');
  events(iEvent).dt_center = [-0.2 0.2 -0.1 0.1]-0.3;
  events(iEvent).v_norm = 35;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [0.985, -0.141, 0.097; 0.152, 0.982, -0.109; -0.080, 0.122, 0.989];    
end
if 0 % Tang 2022
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Tang et al. (2022)';
  events(iEvent).tint = irf.tint('2018-08-27T12:15:36.00Z/2018-08-27T12:15:50.00Z');
  events(iEvent).t_center = irf_time('2018-08-27T12:15:43.50Z','utc>EpochTT');
  events(iEvent).dt_center = [0 0 0.5 -0.5]-0.4;
  events(iEvent).v_norm = 40;
  events(iEvent).sc_id = [1 2 3];
  events(iEvent).lmn = [0.99, -0.12, -0.11; 0.15, 0.96, 0.25; 0.08, -0.26, 0.96];    
end
if 1 % Ergun 2022
  iEvent = iEvent + 1;
  events(iEvent).paper = 'Ergun et al. (2022)';
  events(iEvent).tint = irf.tint('2018-08-27T12:15:36.00Z/2018-08-27T12:15:50.00Z');
  events(iEvent).t_center = irf_time('2018-08-27T12:15:43.50Z','utc>EpochTT');
  events(iEvent).dt_center = [0 0 0.5 -0.5]-0.4;
  events(iEvent).v_norm = 25;
  events(iEvent).sc_id = [1 2 3 4];
  events(iEvent).lmn = [0.910, -0.385,-0.155; 0.415, 0.848, 0.331; 0.004, -0.365, 0.931];
end

nEvents = numel(events);

% Load data
for iEvent = 1%:nEvents
  tint = events(iEvent).tint;
  ic = events(iEvent).sc_id;
  gseB = cell(1,4);
  gseVe = cell(1,4);
  ne = cell(1,4);
  c_eval('gseB{?} = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  c_eval('gseVe{?} = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic)  
  c_eval('ne{?} = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic)

  c_eval('gseR? = gseR?.resample(gseB{?});',1:4)
  c_eval('gseB? = gseB{?};',1:4)
  [Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?','gseB?');
  Jcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data);
  Jcurl.data = Jcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
  Jcurl.time = EpochTT(Jcurl.time); Jcurl.name = '4sc current density';

  events(iEvent).gseB = gseB;
  events(iEvent).gseVe = gseVe;
  events(iEvent).ne = ne;
  events(iEvent).gseJ = Jcurl;
  events(iEvent).gseBav = Bbrst;
end
disp('Done loading data.')

%% Plot data
units = irf_units;
colors = mms_colors('1234');
fontsize = 14;
fontsize_text = 12;

xlim = 49.99*[-1 1];

h = setup_subplots(1,2);
isub = 1;

if 1 % BL
  hca = h(isub); isub = isub + 1;
  for iEvent = 1:nEvents
    event = events(iEvent);
    dy = (1-iEvent)*1; % nT
    t_center = event.t_center;
    dt = event.dt_center;
    dt = dt;
    vnorm = event.v_norm;
    lmn = event.lmn;
    
  
    plot(hca,xlim,[0 0]-dy,'Color',[0.5 0.5 0.5])
    hca.NextPlot = 'add';
    for ic = event.sc_id
      gseB = event.gseB{ic};
      lmnB = gseB*lmn';
      BL0 = lmnB.x.data(1);
      gseVe = event.gseB{ic};
      lmnVe = gseVe*lmn';
  
      n0 = event.ne{ic}.data(1)*1e6; % m^-3
      wpe = sqrt(n0*units.e^2/units.me/units.eps0); % rad/s
      de = units.c/wpe*1e-3; % km
      %disp(sprintf('iEvent = %g, iSc = %g, de = %.0f km',iEvent,ic,de))    
      time = gseB.time;
      time_centered_shifted = time-t_center-dt(ic);
      x = time_centered_shifted*vnorm/de;  
      plot(hca,x,(lmnB.x.data)/BL0-dy,'color',colors(ic,:))
    end
  
    %ht = text(x(1),-dy+1.2,event.paper,'horizontalalignment','right','BackgroundColor',[1 1 0]);
    ht = text(hca,-7,-dy+1.1,event.paper,'horizontalalignment','right','verticalalignment','bottom','fontsize',fontsize_text);
  end
  plot(hca,[0 0],hca.YLim,'Color',[0.5 0.5 0.5])
  hca.NextPlot = 'replaceall';
  hca.XLabel.String = 'N (d_e)';
  hca.YLabel.String = 'B_L/B_{L}^\infty';
  hca.YLabel.String = '$\frac{B_L}{B_L^\infty}$';
  hca.YLabel.FontSize = fontsize + 8;
  hca.YLabel.Interpreter = 'latex';
  
  hca.XLim = xlim;
  %hca.XLim = [-2 2];
end

if 0 % BM
  hca = h(isub); isub = isub + 1;
  for iEvent = 1:nEvents
    event = events(iEvent);
    dy = (1-iEvent)*1; % nT
    t_center = event.t_center;
    dt = event.dt_center;
    %dt = [0.500         -0.05   -0.45   -0.92]+0.1;
    vnorm = event.v_norm;
    lmn = event.lmn;
    
  
    plot(hca,xlim,[0 0]-dy,'Color',[0.5 0.5 0.5])
    hca.NextPlot = 'add';
    for ic = event.sc_id
      gseB = event.gseB{ic};
      lmnB = gseB*lmn';
      BL0 = lmnB.x.data(1);   
      BM0 = lmnB.y.data(1);    
      n0 = event.ne{ic}.data(1)*1e6; % m^-3
      wpe = sqrt(n0*units.e^2/units.me/units.eps0); % rad/s
      de = units.c/wpe*1e-3; % km
      %disp(sprintf('iEvent = %g, iSc = %g, de = %.0f km',iEvent,ic,de))    
      time = gseB.time;
      time_centered_shifted = time-t_center-dt(ic);
      x = time_centered_shifted*vnorm/de;  
      plot(hca,x,(lmnB.y.data)/BL0-BM0/BL0-dy,'color',colors(ic,:))
    end
  
    %ht = text(x(1),-dy+1.2,event.paper,'horizontalalignment','right','BackgroundColor',[1 1 0]);
    %ht = text(-20,-dy+1.2,event.paper,'horizontalalignment','right');
  end
  plot(hca,[0 0],hca.YLim,'Color',[0.5 0.5 0.5])
  hca.NextPlot = 'replaceall';
  hca.XLabel.String = 'N (d_e)';
  hca.YLabel.String = '(B_M-B_M^\infty)/B_{L}^\infty';
  hca.YLabel.String = '$\frac{B_M-B_M^\infty}{|B_L^\infty|}$';
  hca.YLabel.Interpreter = 'latex';
  
  hca.XLim = xlim;
  %hca.XLim = [-2 2];
end

if 0 % veM
  hca = h(isub); isub = isub + 1;
  for iEvent = 1:nEvents
    event = events(iEvent);
    dy = (1-iEvent)*1; % nT
    t_center = event.t_center;
    dt = event.dt_center;
    %dt = [0.500         -0.05   -0.45   -0.92]+0.1;
    vnorm = event.v_norm;
    lmn = event.lmn;
    nSmooth = round(event.gseVe{1}.length/100);
    
  
    plot(hca,xlim,[0 0]-dy,'Color',[0.5 0.5 0.5])
    hca.NextPlot = 'add';
    for ic = event.sc_id
      gseB = event.gseB{ic};
      lmnB = gseB*lmn';
      BL0 = lmnB.x.data(1);
      gseVe = event.gseVe{ic}.smooth(nSmooth);
      gseVeSmooth = event.gseVe{ic}.smooth(nSmooth);
      lmnVe = gseVe*lmn';
      lmnVeSmooth = gseVeSmooth*lmn';
      VeM0 = max(lmnVe.y.data(1));
      VeMax = max(lmnVeSmooth.abs.data);
      lmnVeN = lmnVe*event.ne{ic}.data(1)*1e6;
  
      n0 = event.ne{ic}.data(1)*1e6; % m^-3
      wpe = sqrt(n0*units.e^2/units.me/units.eps0); % rad/s
      de = units.c/wpe*1e-3; % km
      disp(sprintf('iEvent = %g, iSc = %g, de = %.0f km, n_smooth = %g',iEvent,ic,de,nSmooth))    
      time = gseVe.time;
      time_centered_shifted = time-t_center-dt(ic);
      x = time_centered_shifted*vnorm/de;  
      plot(hca,x,sign(BL0)*(lmnVe.y.data-VeM0)/VeMax-dy,'color',colors(ic,:))
    end
  
    %ht = text(x(1),-dy+1.2,event.paper,'horizontalalignment','right','BackgroundColor',[1 1 0]);
    %ht = text(hca,-20,-dy+0.2,event.paper,'horizontalalignment','right');
  end
  plot(hca,[0 0],[-2 nEvents+2],'Color',[0.5 0.5 0.5])
  hca.NextPlot = 'replaceall';
  hca.XLabel.String = 'N (d_e)';
  hca.YLabel.String = 'v_{eM}/V_{eM}^\infty';
  hca.XLim = xlim;
  %hca.XLim = [-2 2];
end

if 1 % jeM
  hca = h(isub); isub = isub + 1;
  for iEvent = 1:nEvents
    event = events(iEvent);
    dy = (1-iEvent)*1; % nT
    t_center = event.t_center;
    dt = event.dt_center;
    %dt = [0.500         -0.05   -0.45   -0.92]+0.1;
    vnorm = event.v_norm;
    lmn = event.lmn;
    nSmooth = round(event.gseVe{1}.length/100);
    
  
    plot(hca,xlim,[0 0]-dy,'Color',[0.5 0.5 0.5])
    hca.NextPlot = 'add';
    for ic = event.sc_id
      gseB = event.gseB{ic};
      lmnB = gseB*lmn';
      BL0 = lmnB.x.data(1);
      gseVe = event.gseVe{ic}.smooth(nSmooth);
      gseVeSmooth = event.gseVe{ic}.smooth(nSmooth);
      lmnVe = gseVe*lmn';
      lmnVeSmooth = gseVeSmooth*lmn';
      VeM0 = max(lmnVe.y.data(1));
      VeMax = max(lmnVeSmooth.abs.data);
      lmnVeN = -lmnVe*event.ne{ic}.data(1)*1e6*units.e*1e3*1e9;
      lmnVeN_smooth = -lmnVe*event.ne{ic}.data(1)*1e6*units.e*1e3*1e9;
      JeM0 = max(lmnVeN.y.data(1));
      JeMax = max(lmnVeN_smooth.abs.data);

      n0 = event.ne{ic}.data(1)*1e6; % m^-3
      wpe = sqrt(n0*units.e^2/units.me/units.eps0); % rad/s
      de = units.c/wpe*1e-3; % km
      disp(sprintf('iEvent = %g, iSc = %g, de = %.0f km, n_smooth = %g',iEvent,ic,de,nSmooth))    
      time = gseVe.time;
      time_centered_shifted = time-t_center-dt(ic);
      x = time_centered_shifted*vnorm/de;  
      plot(hca,x,-sign(BL0)*(lmnVeN.y.data-JeM0)/JeMax-dy,'color',colors(ic,:))
    end
  
    %ht = text(x(1),-dy+1.2,event.paper,'horizontalalignment','right','BackgroundColor',[1 1 0]);
    %ht = text(hca,-20,-dy+0.2,event.paper,'horizontalalignment','right');
  end
  plot(hca,[0 0],[-2 nEvents+2],'Color',[0.5 0.5 0.5])
  hca.NextPlot = 'replaceall';
  hca.XLabel.String = 'N (d_e)';
  hca.YLabel.String = 'B_L^\infty/|B_L^\infty|J_{eM}/J_{eM}^{max}';
  hca.YLabel.String = '$\frac{B_L^\infty}{|B_L^\infty|}\frac{J_{eM}}{\max(J_{eM})}$';
  hca.YLabel.Interpreter = 'latex';
  hca.XLim = xlim;
  %hca.XLim = [-2 2];
end

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).YLabel.FontSize = fontsize + 8;',1:numel(h))
c_eval('h(?).XTick = -100:10:100;',1:numel(h))
c_eval('h(?).XMinorTick = ''on'';',1:numel(h))
c_eval('h(?).XMinorGrid = ''off'';',1:numel(h))
%c_eval('xticklabels = h(?).XTickLabels;',1:numel(h))
c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))

hl = findobj(gcf,'type','line'); c_eval('hl(?).LineWidth = 1;',1:numel(hl))
linkprop(h,{'XLim','YLim'})
h(1).YLim = [-1.5 6.5];


