%% Load DB table
localuser = 'cecilianorgren';
%file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2022.nc'];

file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.nc'];
%data = load(file);
file_csv = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.csv'];
tlim = readtable(file_csv);
tstart = tlim.Var1; tstart.TimeZone = "UTCLeapSeconds";
tstop = tlim.Var2;  tstop.TimeZone = "UTCLeapSeconds";

tstart_ttns = convertTo(tstart, 'tt2000');
tstop_ttns = convertTo(tstop, 'tt2000');
t_duration_s = double((tstop_ttns - tstart_ttns))*1e-9;

%ncdisp(file)

info = ncinfo(file);
vars = {info.Variables.Name};
nvars = numel(vars);

clear db
for ivar = 1:nvars
  db.(vars{ivar}) = ncread(file,vars{ivar});
end

db_table_ff = struct2table(db);
% Add stop time of slow
db_table_ff = addvars(db_table_ff,tstart_ttns,tstop_ttns,t_duration_s,'After','time');


t0_ff = EpochTT('2017-05-05T19:41:44.790324');
% time: microseconds since 2017-05-05T19:41:44.790324
time_ff_ttns = t0_ff.ttns + db_table_ff.time*1e3;
db_table_ff.time = time_ff_ttns; % rewrite time in ttns


t0_df = EpochTT('2017-05-19T03:06:44.458185978');
% t_df: nanoseconds since  2017-05-19T03:06:44.458185978
time_df_ttns = double(t0_df.ttns + int64(db_table_ff.t_df));
time_df_ttns(~db_table_ff.is_df) = NaN;
db_table_ff.t_df = time_df_ttns; % rewrite time in ttns


%% Loop events and plot brief overview
units = irf_units;

mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');

tint = irf.tint('2017-07-09T17:30:50.000Z/2017-07-09T17:33:50.000Z');

iDF = 0;
doPrint = 0*1;
doPlot = 0;
vecK = 1:5;

% Turn off irf._log warning messages
irf.log('off')
% Parameters
ic = 1;
nMovMean = 1; % Number of distributions for the moving average. Use the same for finding the energy limt and the gmm
mask_nMovMean = [7 7]; mask_Counts = 5;
nMP = 10000; % Number of macroparticles
K = 1:6;  
reduction.metric = 'l1'; % for reduction
reduction.threshold = [0.1 0.2 0.3 0.4]; % for reduction
reduction2.metric = 'l1'; % for reduction
reduction2.threshold = 0.02:0.02:0.3; % for reduction

doPrint = true;
  str_event = sprintf('Ev%03.0f_nMovMean=%g_mask=%.0f_%.0f',iDF,nMovMean,mask_nMovMean(1),mask_nMovMean(1));
  
  
  disp('------------------------------------------------------------')
  disp("iDF = " + iDF + ...
    ", t1 = " + tint(1).utc('yyyy-mm-ddTHH:MM:SS') + ...
    ", t2 = " + tint(2).utc('yyyy-mm-ddTHH:MM:SS'))
  
  
  disp('>> Loading data.')
  try
    c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
    c_eval('gseVi = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
    c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);

    c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
    c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
    c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
    iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;
  catch
     warning('Error loading data for ic=%g. Continuing to next event.',ic)
     %continue
  end
  %% Prep data
  disp('>> Preparing data.')
  PD_use = iPDist; % iPDist
  %mask_nMovMean = [7 7]; mask_Counts = 5;
  %run_metadata.PD.mask_nMovMean = mask_nMovMean;
  %run_metadata.PD.mask_nCounts = mask_Counts;

  mat = iPDist_counts.find_low_counts('counts',mask_Counts,'nMovMean',mask_nMovMean,'output','mat');
  PD = PD_use.mask('energy','mat',mat);  
  PD = PD.movmean(nMovMean); % see how this affects the results
  PD = PD(1:nMovMean:PD.length); % This step is done just to reduce the computations for now, need to check it it is needded for noise
  PDc = iPDist_counts.movmean(nMovMean);
  PDc = PDc(1:nMovMean:iPDist_counts.length);

  MP = PD.macroparticles('ntot',nMP,'skipzero',1,'scpot',scPot.resample(PD),'v_init','v2','MParts');

  fz = PD.reduce('1D',[0 0 1]);



  %% Run GMM
  disp('>> Running fit.')
  options = statset('MaxIter',150,'TolFun',1e-5);
  GM = GMDist.fit(MP,K,'start','previous','Options',options,'SharedCovariance',false,'RegularizationValue',0.1);

  %% Reduce mixture
  disp('>> Reducing mixture.')
  nCompPreMerge = GM.nComp;
  nCompMax = max(nCompPreMerge(:));
  nCompMin = min(nCompPreMerge(:));
  mGM = cell(1,numel(reduction.threshold));
  
  fprintf(sprintf('Caching merge error, %s.\n',reduction.metric))  
  GM = GM.cacheMergeErr('metric',reduction.metric);        % the slow step, once
  fprintf(sprintf('Applying thresholds.\n',reduction.metric))
  for ithr = 1:numel(reduction.threshold)
    fprintf(sprintf('%s<%.3f\n',reduction.metric,reduction.threshold(ithr)))
    mGM{ithr} = GM.merge('threshold',reduction.threshold(ithr),'metric',reduction.metric);
  end

  for ithr = 1:numel(reduction2.threshold)
    fprintf(sprintf('%s<%.3f\n',reduction2.metric,reduction2.threshold(ithr)))
    mGM2{ithr} = GM.merge('threshold',reduction2.threshold(ithr),'metric',reduction2.metric);
  end

  alpha = 1.5; floor = 0.05;
  thr6 = GM.slice(':',2).weightThreshold('alpha',alpha,'floor',floor);        % [nT x 1]
  mGM3 = GM.merge('threshold',thr6,'metric','l1');  
  mGM3.userdata.threshold_alpha = alpha;
  mGM3.userdata.threshold_floor = floor;
  tsThr = irf.ts_scalar(GM.time,thr6);

  % tic
  % for ithr = 1:numel(reduction.threshold)
  %   disp(sprintf('%s<%.3f',reduction.metric,reduction.threshold(ithr)))
  %   % One could probably calculate the merge metric once and then apply
  %   % different thresholds to that, not sure how much time it would save.
  %   % Yes, the bestPartition is much slower than merge, factor 30 or so.
  % 
  %   % Ideally one should have the option to combine these two, and bestM, 
  %   % bestP, bestErr would be saved in the object.
  %   tic
  % 
  %   [bestM, bestP, bestErr] = GM.bestPartition(reduction.threshold(ithr),'metric',reduction.metric);
  %   toc
  %   % Here, threshold is only to pass the info of the parameter, not elegant 
  %   % and unstable because one could pass the wrong value
  %   tic
  %   mGM{ithr} = GM.merge(bestP,'threshold',reduction.threshold(ithr)); 
  %   toc
  % end
  % toc

  %% Figure, OMNIDEFlux, f(vz), BIC, K, M(threshold)
  figure(106)
  
  K = 2;
  sortBy = 'T';
  G = GM.slice(8:GM.nT,K).sortBy(sortBy);
  PD_ = PD(8:GM.nT);

  str_print = sprintf('%s_M_dep_on_threshold',str_event);  
  h = irf_plot(5);

  hca = irf_panel('B');
  irf_plot(hca,gseB1)
  hca.ColorOrder = mms_colors('xyz');
  hca.YLabel.String = 'B (nT)';

  hca = irf_panel('omni fpi');
  [hs,hcb] = irf_spectrogram(hca,PD.deflux.omni.specrec);

  if 0
  hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,fz.specrec)
  hca.CLim = prctile(hca.Children.CData(:),[1 99]);
  end

  %hca = irf_panel('BIC'); % it is easier to see counter-streaming populations in this one
  %irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec('p_label','BIC/max(BIC)'),'lin')
  
  hca = irf_panel('T'); % it is easier to see counter-streaming populations in this one
  irf_plot(hca,G.tsT('type','iso'),'comp')
  hca.YScale = 'log';
  if 1
  hold(hca,'on')
  %[Tf,dE,iBin] = GM.slice(':',6).floorTBin(PD);
  tsTFloor = G.tsFloorTBin(PD_);
  irf_plot(hca,tsTFloor,'comp','linestyle','--')
  hold(hca,'off')
  end
  if strcmp(hca.YScale,'log'); hca.YTick = 10.^(-10:10); end
  irf_zoom(hca,'y')
  hca.YLabel.String = 'T (eV)';

  hca = irf_panel('N'); % it is easier to see counter-streaming populations in this one
  irf_plot(hca,G.tsN,'comp')
  %hca.YScale = 'log';
  hca.YLabel.String = 'n (cm^{-3})';

  hca = irf_panel('v'); % it is easier to see counter-streaming populations in this one
  ts = G.tsV; ts = cellfun(@(x){x.z},ts);
  irf_plot(hca,ts,'comp')
  %hca.YScale = 'log';
  hca.YLabel.String = 'v_z (km/s)';
  

  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h(?).Layer = ''top'';',1:numel(h))
  irf_plot_axis_align
  irf_zoom(h,'x',G.time)
  irf_timeaxis(h)
  h(end).XTickLabelRotation = 0;
  hl = findobj(gcf,'type','line'); c_eval('hl(?).LineWidth = 1.;',1:numel(hl))

  c_eval('h(?).FontSize = 10;',1:numel(h))
  irf_plot_axis_align(h)

  if doPrint, cn.print(str_print); end

  

  %% Figure, Bin heating
  figure(17007)
  
  K = 2;
  sortBy = 'T';
  G = GM.slice(1:GM.nT,K).sortBy(sortBy);
  PD_ = PD(1:GM.nT);

  str_print = sprintf('%s_M_dep_on_threshold',str_event);  
  h = irf_plot(2);
  c_eval('h(?).Position(2)=h(?).Position(2)+0.04;',1:2)

  if 0
  hca = irf_panel('B');
  irf_plot(hca,gseB1)
  hca.ColorOrder = mms_colors('xyz');
  hca.YLabel.String = 'B (nT)';
  end

  hca = irf_panel('omni fpi');
  [hs,hcb] = irf_spectrogram(hca,PD.deflux.omni.specrec);

  if 0
  hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,fz.specrec)
  hca.CLim = prctile(hca.Children.CData(:),[1 99]);
  end

  %hca = irf_panel('BIC'); % it is easier to see counter-streaming populations in this one
  %irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec('p_label','BIC/max(BIC)'),'lin')
  
  hca = irf_panel('T'); % it is easier to see counter-streaming populations in this one
  ts = G.tsT('type','iso');
  irf_plot(hca,ts{1},'comp')
  hca.ColorOrder = mms_colors('21');
  hca.YScale = 'log';
  if 1
  hold(hca,'on')
  %[Tf,dE,iBin] = GM.slice(':',6).floorTBin(PD);
  tsTFloorIso = G.sortBy('T').tsFloorTBin(PD_,'type','iso');
  tsTFloorRad = G.sortBy('T').tsFloorTBin(PD_,'type','radial');
  irf_plot(hca,tsTFloorIso{1},'comp','linestyle','--')
  %irf_plot(hca,tsTFloorRad{1},'comp','linestyle','-')
  hold(hca,'off')
  end
  if strcmp(hca.YScale,'log'); hca.YTick = 10.^(-10:10); end
  irf_zoom(hca,'y')
  hca.YLabel.String = 'T (eV)';
  irf_legend(hca,{'T_{coldest}','T_{floor}'},[0.98 0.1])

  
  h(1).CLim = [3.8 6.5];

  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h(?).Layer = ''top'';',1:numel(h))
  irf_plot_axis_align
  irf_zoom(h,'x',G.time)
  irf_timeaxis(h)
  h(end).XTickLabelRotation = 0;
  hl = findobj(gcf,'type','line'); c_eval('hl(?).LineWidth = 1.;',1:numel(hl))

  c_eval('h(?).FontSize = 10;',1:numel(h))
  irf_plot_axis_align(h)

  if doPrint, cn.print(str_print); end

  

 
  %% Figure, OMNIDEFlux, f(vz), BIC, K, M(threshold)
  figure(101)
  str_print = sprintf('%s_M_dep_on_threshold',str_event);
  h = irf_plot(5);

  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)
  hca.CLim = [3.8 6.5];

  hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,fz.specrec)
  %hca.CLim = prctile(hca.Children.CData(:),[1 99]);
  
  hca = irf_panel('BIC'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec('p_label','BIC/max(BIC)'),'lin')
  
  for K = 1:2
    hca = irf_panel(sprintf('GM fz, K=%g',K)); % it is easier to see counter-streaming populations in this one
    irf_spectrogram(hca,GM.slice(':',K).tsReduced(3).specrec,'log')
    %hca.CLim = prctile(hca.Children.CData(:),[1 99]);
    %hca.CLim = prctile(hca.Children.CData(:),[1 99]);
  end

  hlinks = linkprop(h([2 4 5]),{'XLim','CLim'});

  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h(?).Layer = ''top'';',1:numel(h))
  irf_plot_axis_align
  irf_zoom(h,'x',PD.time)
  irf_timeaxis(h)
  h(end).XTickLabelRotation = 0;

  if doPrint, cn.print(str_print); end

  %% Figure, M(K,thr)
  figure(103)  
  
  str_print = char(str_event + "_Scatter_B_nComp_vs_thr_and_K");

  oG = GM.slice(':',oK);
  edgesM = 0.5:(max(K)+1);
  centerM = K;
  NM = nan(numel(reduction.threshold),numel(K),numel(edgesM)-1);
  for ithr = 1:numel(reduction.threshold)
    for iK = 1:numel(K)
      nComp = mGM{ithr}.slice(':',K(iK)).nComp;
      NM(ithr,iK,:) = histcounts(nComp,edgesM);
    end
  end
  
  [h1,h] = initialize_combined_plot('topbottom',1,2,2,0.3,'vertical');
   
  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)     


  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h1(?).Layer = ''top'';',1:numel(h1))
  irf_plot_axis_align
  irf_zoom(h1,'x',PD.time)
  irf_timeaxis(h1)
  h1(end).XTickLabelRotation = 0;
  
  
    

  isub = 1;
  for ithr = 1:numel(reduction.threshold)
    str_title = upper(mGM{ithr}.misc.metric) + " < " + mGM{ithr}.misc.threshold;
    hca = h(isub); isub = isub + 1;
    plotd = squeeze(NM(ithr,:,:));
    bar(hca,centerM,plotd',1)
    hca.XLabel.String = 'M';
    hca.YLabel.String = 'counts';
    hca.Title.String = str_title;
    irf_legend(hca,("K = " + K)',[0.98 0.98])
     
    compact_panels(h,0.07,0.08,1)
    c_eval('h(?).Box = ''on'';',1:numel(h))
    c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
    c_eval('h(?).ColorOrder = mms_colors(''matlab'');',1:numel(h))
  
      
  end
  if doPrint, cn.print(str_print); end
   
  %% Figure, scatterplot of T,T; n,T
  figure(107)
  K = 2;
  sortBy = 'T';
  G = GM.slice(8:GM.nT,K);
  params = G.sortBy(sortBy).byNComp();

  h = setup_subplots(1,2);  
  isub = 1;


  hca = h(isub); isub = isub + 1;
  scatter(hca,params.n,params.T,12,'filled');
  hca.XLabel.String = 'n (cm^{-3})';
  hca.YLabel.String = 'T (eV)';
  %hca.Title.String = sprintf('K = %d',K);
  grid(hca,'on');
  %addBackground(hca,@(x,y)x./y,[0.001 0.1 1 10 100],'ShowText',1)
  %hca.XScale = 'log';
  %hca.YScale = 'log';

  hca = h(isub); isub = isub + 1;
  scatter(hca,params.T(:,1),params.T(:,2),12,'filled');
  hca.XLabel.String = 'T_{cold} (eV)';
  hca.YLabel.String = 'T_{hot} (eV)';
  %hca.Title.String = sprintf('K = %d',K);
  grid(hca,'on');
  hca.XScale = 'log';
  hca.YScale = 'log';
  addBackground(hca,@(x,y)x./y,[0.01 0.1 1 10 100],'ShowText',1)
  
  hca.XTick = 10.^(0:5);
  %hca.YTick = 10.^(0:5);


  if 0
  hca = h(isub); isub = isub + 1;
  histNorm = 'count';
  histogram(hca,params.T(:,1),logspace(1,5,25),'Normalization',histNorm)
  hold(hca,'on')
  histogram(hca,params.T(:,2),logspace(1,5,25),'Normalization',histNorm)
  hold(hca,'off')
  hca.XScale = 'log';
  hca.XTick = 10.^(0:5);
  end

  c_eval('h(?).Box = ''on'';',1:numel(h))
  

  %% Figure, M(K,thr)
  figure(103)  
  
  str_print = char(str_event + "_Scatter_B_nComp_vs_thr_and_K");

  oG = GM.slice(':',oK);
  edgesM = 0.5:(max(K)+1);
  centerM = K;
  NM = nan(numel(reduction.threshold),numel(K),numel(edgesM)-1);
  for ithr = 1:numel(reduction.threshold)
    for iK = 1:numel(K)
      nComp = mGM{ithr}.slice(':',K(iK)).nComp;
      NM(ithr,iK,:) = histcounts(nComp,edgesM);
    end
  end
  
  [h1,h] = initialize_combined_plot('topbottom',1,2,2,0.3,'vertical');
   
  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)     


  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h1(?).Layer = ''top'';',1:numel(h1))
  irf_plot_axis_align
  irf_zoom(h1,'x',PD.time)
  irf_timeaxis(h1)
  h1(end).XTickLabelRotation = 0;
  
  
    

  isub = 1;
  for ithr = 1:numel(reduction.threshold)
    str_title = upper(mGM{ithr}.misc.metric) + " < " + mGM{ithr}.misc.threshold;
    hca = h(isub); isub = isub + 1;
    plotd = squeeze(NM(ithr,:,:));
    bar(hca,centerM,plotd',1)
    hca.XLabel.String = 'M';
    hca.YLabel.String = 'counts';
    hca.Title.String = str_title;
    irf_legend(hca,("K = " + K)',[0.98 0.98])
     
    compact_panels(h,0.07,0.08,1)
    c_eval('h(?).Box = ''on'';',1:numel(h))
    c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
    c_eval('h(?).ColorOrder = mms_colors(''matlab'');',1:numel(h))
  
      
  end
  if doPrint, cn.print(str_print); end
   
  %% Figure, M(K,thr) 2
  figure(104)  
  
  str_print = char(str_event + "_HistogramStacked__M_vs_thr_and_K");

  oG = GM.slice(':',oK);
  nCompAll = cellfun(@(x)x.nComp,mGM2,'UniformOutput',false); 
  uM = unique(cat(1,nCompAll{:}));
  %edgesErr = 0:0.02:0.3;
  %centerErr = edgesErr(1:end-1)+0.5*diff(edgesErr);
  histNorm = 'pdf';
  NM = nan(numel(reduction2.threshold),numel(K),numel(K));
  for ithr = 1:numel(reduction2.threshold)
    for iK = 1:numel(K)
      nComp = mGM2{ithr}.slice(':',K(iK)).nComp;
      err = mGM2{ithr}.misc.merged_partition.err(:,iK);      
      for iM = 1:numel(K(1:iK))        
        %N_ = histcounts(err(nComp==K(iM)),edgesErr,'Normalization',histNorm);
        NM(ithr,iM,iK) = numel(find(nComp==K(iM)));
      end
    end
  end
  
  [h1,h] = initialize_combined_plot('topbottom',1,2,3,0.3,'vertical');
   
  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)     


  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h1(?).Layer = ''top'';',1:numel(h1))
  irf_plot_axis_align
  irf_zoom(h1,'x',PD.time)
  irf_timeaxis(h1)
  h1(end).XTickLabelRotation = 0;
      

  isub = 1;
  for iK = 1:numel(K)
    str_title = "K=" + K(iK);
    hca = h(isub); isub = isub + 1;
    plotd = squeeze(NM(:,:,iK));
    %bar(hca,centerErr,plotd',1)
    bar(hca,reduction2.threshold,plotd'/GM.nT,1,'stacked','FaceAlpha',0.5)
    % for iM = 1:numel(uM)
    %   bar(hca,reduction2.threshold,plotd(:,iM),1,'FaceAlpha',0.5)
    %   if iM==1, hold(hca,'on'); end
    % end
    hold(hca,'off');
    hca.YLim = [0 1];
    hca.XLabel.String = upper(reduction2.metric);
    hca.YLabel.String = histNorm;
    hca.Title.String = str_title;
    irf_legend(hca,("M = " + K(1:iK)'),[1.02 0.98])
     
    compact_panels(h,0.07,0.08,1)
    c_eval('h(?).Box = ''on'';',1:numel(h))
    c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
    c_eval('h(?).ColorOrder = mms_colors(''matlab'');',1:numel(h))
  
      
  end
  if doPrint, cn.print(str_print); end
 
  %% Figure, density and temperature for different thresholds and M=2,3
  figure(102)
  sortBy = 'T';  
  oK = 6;
  oG = GM.slice(':',oK);

  for ithr = 1%1:numel(reduction.threshold)
    G = mGM{ithr}.slice(':',oK);
    uM = unique(G.nComp);
    params = G.sortBy(sortBy).byNComp();
    dims = {'x','y','z'};
     
    
    [h1,h] = initialize_combined_plot('topbottom',1,2,4,0.3,'vertical');
  
    hca = irf_panel('M');
    irf_plot(hca,{oG.tsNComp,G.tsNComp},'comp')
    hca.YLabel.String = 'nComp';
    hca.YLim = [0 oK+1];
    irf_legend(hca,{'K',sprintf('M (%s<%.2f)',upper(mGM{ithr}.misc.metric),mGM{ithr}.misc.threshold)},[0.98 0.10])
    
    irf_zoom(h1,'x',PD.time)
    irf_timeaxis(h1)
    h1(end).XTickLabelRotation = 0;
  
    for M = 3%uM(:)'    
      iM = find(uM==M);    
      nData = numel(params(iM).entry);
      str_var = "K=" + oK + "_" + upper(mGM{ithr}.misc.metric) + "<" + mGM{ithr}.misc.threshold+ "_M=" + M  + "_sortBy=" + sortBy;
      str_print = char(str_event + "_Scatter_A_" + str_var);
      str_title = "K=" + oK + ", M=" + M + " (" + upper(mGM{ithr}.misc.metric) + "<" + mGM{ithr}.misc.threshold + "), sortBy=" + sortBy + ", nData = " + nData;
    
  
      isub = 1;
    
      hca = h(isub); isub = isub + 1;
      scatter(hca,params(iM).n,params(iM).T)
      hca.XLabel.String = 'n (cc)';
      hca.YLabel.String = 'T (eV)';
      %hca.Title.String = str_title;
    
      hca = h(isub); isub = isub + 1;
      scatter(hca,params(iM).n(:,end),params(iM).n(:,1:end-1))
      hca.XLabel.String = 'n_{last} (cc)';
      hca.YLabel.String = 'n_{others} (cc)';
      addBackground(hca,@(x,y)x./y,[1:10],'ShowText',1)
      
      hca = h(isub); isub = isub + 1;
      scatter(hca,params(iM).T(:,end),params(iM).T(:,1:end-1))
      hca.XLabel.String = 'T_{last} (eV)';
      hca.YLabel.String = 'T_{others} (eV)';
      %addBackground(hca,@(x,y)x./y,[1 2 5 10 20 50 100 200 500 1000])
      addBackground(hca,@(x,y)x./y,[1 10 100],'ShowText',1)
      hca.XScale = 'log';
      hca.YScale = 'log';
      
      
      hca = h(isub); isub = isub + 1;
      edgesT = linspace(1,5,25);
      centerT = edgesT(1:end-1) + 0.5*diff(edgesT);
      normHist = 'pdf';
      fun = @(x)histcounts(log10(x)',edgesT,'normalization',normHist);
      pdf = rowfun(fun,table(params(iM).T')); pdf = pdf.Var1;
      for irow = 1:size(pdf,1)
        hcbar = bar(hca,centerT,pdf(irow,:),1,'facealpha',0.5);
        if irow == 1; hold(hca,'on'), end
      end  
      hold(hca,'off')
      hca.ColorOrder = mms_colors('matlab');
      hca.XLabel.String = 'log_{10} T (eV)';
      hca.XTick = -1:6;
      hca.XTickLabel = "10^" + hca.XTick;
      hca.YLabel.String = normHist;
    
    
      hca = h(isub); isub = isub + 1;
      xplot = params(iM).T;
      yplot = params(iM).Uk;
      scatter(hca,xplot,yplot)      
      hca.YLabel.String = 'U_k (eV)';
      hca.XLabel.String = 'T (eV)';
      addBackground(hca,@(x,y)x./y,10.^(-2:2),'ShowText',1)
      %hca.XScale = 'log';      
      %hca.YScale = 'log';
      hca.XLim = prctile(xplot(:),[1 99]);
      hca.YLim = prctile(yplot(:),[1 99]);
      
    
      for idim = 1:3
        hca = h(isub); isub = isub + 1;
        xplot = params(iM).V(:,:,idim);
        yplot = log10(params(iM).T);
        scatter(hca,xplot,yplot)
        hca.YLabel.String = 'T (eV)';
        hca.XLabel.String = sprintf('v_%s (km/s)',dims{idim});      
        hca.XLim = prctile(xplot(:),[1 99]);
        hca.YLim = prctile(yplot(:),[1 99]);
      end

      h1.Title.String = str_title;
      h1.Title.Position(2)=-4;
      
      hs = findobj(gcf,'type','scatter');
      %c_eval('hs(?).Marker = ''o''; hs(?).MarkerFaceColor = "flat";',1:numel(hs))
      compact_panels(h,0.07,0.05,1)
      c_eval('h(?).Box = ''on'';',1:numel(h))
      c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
      c_eval('h(?).ColorOrder = mms_colors(''matlab'');',1:numel(h))
  
      if doPrint, cn.print(str_print); end
    end
  end
  
  %% Figure, varying threshold
  thrRefK = 6;
  plotK = 6;
  alpha = 1.0; floor = 0.10; % in terms of weight
  thr = GM.slice(':',thrRefK).weightThreshold('alpha',alpha,'floor',floor);        % [nT x 1]
  mGM3 = GM.merge('threshold',thr,'metric','l1');  
  mGM3.userdata.threshold_alpha = alpha;
  mGM3.userdata.threshold_floor = floor;
  tsThr = irf.ts_scalar(GM.time,thr);
  
  [cold,info_cold] = mGM3.cold('maxColdEV',500,'minGapLog10',0.7);
  tsNCold = irf.ts_scalar(mGM3.time,cellfun(@(x) sum(x),cold(:,plotK)));
  tsNHot = irf.ts_scalar(mGM3.time,cellfun(@(x) sum(~x),cold(:,plotK)));
%
  figure(105)
  str_print = sprintf('%s_varying_threshold_alpha=%.2f_floor=%.2f',str_event,alpha,floor);
  h = irf_plot(8);

  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)

  hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,fz.specrec)
  hca.CLim = prctile(hca.Children.CData(:),[1 99]);
    
  hca = irf_panel('nComp before merge');
  [~,hcb] = irf_spectrogram(hca,GM.tsNComp.specrec('f',nCompMin:nCompMax, 'f_label',{'K'},'p_label',{'M'}),'lin');
  hcb.YTick = 0:nCompMax;
  hca.CLim = [0 nCompMax];

  hca = irf_panel('nComp after merge');
  [~,hcb] = irf_spectrogram(hca,mGM3.tsNComp.specrec('f',nCompMin:nCompMax, 'f_label',{'K'},'p_label',{'M'}),'lin');
  hcb.YTick = 0:nCompMax;
  hca.CLim = [0 nCompMax]; 

  if 0 % timevarying threshold
    hca = irf_panel('timevarying threshold');
    irf_plot(hca,tsThr)
    hca.YLabel.String = 'thr';
    irf_legend(hca,{"alpha = " + mGM3.userdata.threshold_alpha + ", floor = " + mGM3.userdata.threshold_floor},[0.98 0.98],'color','k')
  end
  if 1 % weights and timevarying threshold
    hca = irf_panel('ComponentProportion after merge and timevarying threshold');
    hl = irf_plot(hca,mGM3.slice(':',plotK).sortBy('n','direction','descend').tsProp,'comp');
    c_eval('hl.Children(?).Marker = "*"; hl.Children(?).LineStyle = "none";',1:numel(hl.Children))
    hold(hca,'on')
    irf_plot(hca,tsThr,'k-')
    hold(hca,'off')
    hca.YLabel.String = 'weight';     
    irf_legend(hca,cat(2,"M=" + (1:max(mGM3.slice(':',plotK).nComp)),"threshold")',[1.02 0.98])
  end

  if 0 % density after merge'
    hca = irf_panel('density after merge');
    hl = irf_plot(hca,mGM3.slice(':',plotK).sortBy('n','direction','descend').tsN,'comp');
    c_eval('hl.Children(?).Marker = "*"; hl.Children(?).LineStyle = "none";',1:numel(hl.Children))
    hca.YLabel.String = 'n (cc)';
    %hca.YScale = 'log';
    %hca.YLim(1) = 10^-3;
  end
  if 0 % ComponentProportion after merge'
    hca = irf_panel('ComponentProportion after merge');
    hl = irf_plot(hca,mGM3.slice(':',plotK).sortBy('n','direction','descend').tsProp,'comp');
    c_eval('hl.Children(?).Marker = "*"; hl.Children(?).LineStyle = "none";',1:numel(hl.Children))
    hca.YLabel.String = 'weight';
    %hca.YScale = 'log';
    %hca.YLim(1) = 10^-3;
  end


  hca = irf_panel('reduced f after merge');
  irf_spectrogram(hca,mGM3.slice(':',plotK).tsReduced(3).specrec,'log')
  %hca.YLabel.String = 'n (cc)';
  %hca.YScale = 'log';
  %hca.YLim(1) = 10^-3;


  hca = irf_panel('temeprature after merge');
  hl = irf_plot(hca,mGM3.slice(':',plotK).sortBy('n','direction','descend').tsT('type','iso'),'comp');
  c_eval('hl.Children(?).Marker = "*"; hl.Children(?).LineStyle = "none";',1:numel(hl.Children))
  hca.YLabel.String = 'T (eV), sorted';
  hca.YScale = 'log';
  %hca.YLim(1) = 10^-3;
  irf_zoom(hca,'y')
  
  hca = irf_panel('NCold, NHot');
  irf_plot(hca,{tsNCold,tsNHot},'comp')
  hca.YLabel.String = 'Count';
  irf_legend(hca,{'N_{comp,cold}','N_{comp,hot}'},[0.98 0.98])


  hlinks = linkprop(h([2 6]),{'CLim','YLim'});
  h(2).CLim = [-5 -0.5];

  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h(?).Layer = ''top'';',1:numel(h))
  irf_plot_axis_align
  irf_zoom(h,'x',PD.time)
  %irf_zoom(h,'y')
  irf_timeaxis(h)
  h(end).XTickLabelRotation = 0;

  if doPrint, cn.print(str_print); end


irf.log('on')