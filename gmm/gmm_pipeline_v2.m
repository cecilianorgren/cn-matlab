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

% BBF/DF database, to extract times
db_table_df = db_table_ff(db_table_ff.is_df==1,:);
nDF = numel(db_table_df.time);


iDFs = [9 128];
iDFs = 217:nDF;
%iDFs = 87;
iDFs = 17;23:nDF;
iDFs = 101:nDF;
iDFs = [1 6 9 10 13 21 25 35 38 45 46 47 52 55 56 57 58 60 61 68 69 70 71 78 79 81 83 87 91];
iDFs = [92 95 98 101 107 112 115 121 ];
iDFs = 37;
iDFs = 57;
iDFs=9;
iDFs = 10;
iDFs = 16;
doPrint = 1;
doPlot = 0;
vecK = 1:5;

% Turn off irf._log warning messages
irf.log('off')
% Parameters
ic = 1;
nMovMean = 5; % Number of distributions for the moving average. Use the same for finding the energy limt and the gmm
mask_nMovMean = [5 5]; mask_Counts = 4;
nMP = 10000; % Number of macroparticles
K = 1:6;  
reduction.metric = 'l1'; % for reduction
reduction.threshold = [0.1 0.2 0.3 0.4 0.5]; % for reduction
doPrint = true;

% Run through events
for iDF = 7:50%37;iDFs;
  str_event = sprintf('Ev%03.0f',iDF);
  
  % Define time
  t0 = EpochTT(db_table_df.time(iDF));
  T =  db_table_df.t_duration_s(iDF);  
  t_ff_start = t0;
  t_ff_stop = t_ff_start + T;
  tint = [t_ff_start t_ff_stop] + 5*[-1 1];
  %tint = tint + [-90 0];
  tDF = EpochTT(int64(db_table_df.t_df(iDF)));
  
  disp('------------------------------------------------------------')
  disp("iDF = " + iDF + ...
    ", t1 = " + tint(1).utc('yyyy-mm-ddTHH:MM:SS') + ...
    ", t2 = " + tint(2).utc('yyyy-mm-ddTHH:MM:SS'))
  
  
  disp('>> Loading data.')
  try
    c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
    c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
    c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
    iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;
  catch
     warning('Error loading data for ic=%g. Continuing to next event.',ic)
     continue
  end
  %% Prep data
  disp('>> Preparing data.')
  PD_use = iPDist; % iPDist
  %mask_nMovMean = [7 7]; mask_Counts = 5;
  %run_metadata.PD.mask_nMovMean = mask_nMovMean;
  %run_metadata.PD.mask_nCounts = mask_Counts;

  mat = iPDist_counts.find_low_counts('counts',mask_Counts,'nMovMean',mask_nMovMean,'output','mat');
  PD = PD_use.mask('energy','mat',mat);  
  PD = PD.movmean(nMovMean);
  PD = PD(1:nMovMean:PD.length); % This step is done just to reduce the computations for now, need to check it it is needded for noise
  PDc = iPDist_counts.movmean(nMovMean);
  PDc = PDc(1:nMovMean:iPDist_counts.length);

  MP = PD.macroparticles('ntot',nMP,'skipzero',1,'scpot',scPot.resample(PD),'v_init','v2','MParts');

  %fz = PD.reduce('1D',[0 0 1]);

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
  tic
  for ithr = 1:numel(reduction.threshold)
    disp(sprintf('%s<%.3f',reduction.metric,reduction.threshold(ithr)))
    % One could probably calculate the merge metric once and then apply
    % different thresholds to that, not sure how much time it would save.
    % Yes, the bestPartition is much slower than merge, factor 30 or so.

    % Ideally one should have the option to combine these two, and bestM, 
    % bestP, bestErr would be saved in the object.
    tic
    [bestM, bestP, bestErr] = GM.bestPartition(reduction.threshold(ithr),'metric',reduction.metric);
    toc
    % Here, threshold is only to pass the info of the parameter, not elegant 
    % and unstable because one could pass the wrong value
    tic
    mGM{ithr} = GM.merge(bestP,'threshold',reduction.threshold(ithr)); 
    toc
  end
  toc

  %% Figure, merging, how M depends on errL1 threshold
  str_figure = sprintf('%s_M_dep_on_threshold',str_event);
  h = irf_plot(3+numel(mGM));

  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,PD.deflux.omni.specrec)

  hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
  irf_spectrogram(hca,fz.specrec)
  hca.CLim = prctile(hca.Children.CData(:),[1 99])

  hca = irf_panel('nComp before merge');
  [~,hcb] = irf_spectrogram(hca,GM.tsNComp.specrec('f',nCompMin:nCompMax, 'f_label',{'nComp';'pre-merge'},'p_label',{'nComp';'pre-merge'}),'lin');
  hcb.YTick = 0:nCompMax;
  hca.CLim = [0 nCompMax];

  for ii = 1:numel(mGM)
    G = mGM{ii};
    hca = irf_panel(sprintf('nComp after merge ithr=%.2f',threshold(ii)));
    specrec = G.tsNComp.specrec( ...
      'f',nCompMin:nCompMax, ...
      'f_label',{'nComp';'pre-merge'}, ...
      'p_label',{'nComp';'post-merge';sprintf('thr = %.2f',G.misc.threshold)});
    [~,hcb] = irf_spectrogram(hca,specrec,'lin');
    hcb.YTick = 0:nCompMax;
    hca.CLim = [0 nCompMax];
  end


  colormap([flipdim(irf_colormap('Spectral'),1)])
  c_eval('h(?).Layer = ''top'';',1:numel(h))
  irf_plot_axis_align
  irf_zoom('x',PD.time)
  irf_timeaxis(h)
  h(end).XTickLabelRotation = 0;

  if doPrint, cn.print(str_figure); end

end

irf_log on