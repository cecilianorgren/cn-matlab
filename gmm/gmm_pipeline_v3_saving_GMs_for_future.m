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


% BBF/DF database, to extract times
iFFs = 1:size(db_table_ff,1);
iFF_with_DF = find(db_table_ff.is_df==1);
db_table_df = db_table_ff(iFF_with_DF,:);
nDF = numel(db_table_df.time);
db_table = db_table_ff;

%% Loop events and plot brief overview
units = irf_units;

mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');



iDFs = [9 128];
%iDFs = 217:nDF;
%iDFs = 87;
%iDFs = 17;23:nDF;
%iDFs = 101:nDF;
iDFs = [1 6 9 10 13 21 25 35 38 45 46 47 52 55 56 57 58 60 61 68 69 70 71 78 79 81 83 87 91];
iDFs = [92 95 98 101 107 112 115 121 ];
iDFs = 37;
iDFs = 57;
iDFs=9;
iDFs = 10;
iDFs = 16;




% Turn off irf._log warning messages
irf.log('off')

% Save path
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v2/'];

% Parameters
ic = 1;
nMovMean = 5; % Number of distributions for the moving average. Use the same for finding the energy limt and the gmm
mask_nMovMean = [7 7]; mask_Counts = 5;
nMP = 10000; % Number of macroparticles
K = 1:7;
reduction.metric = 'l1'; % for reduction
reduction.threshold = [0.1 0.2 0.3 0.4]; % for reduction
reduction2.metric = 'l1'; % for reduction
reduction2.threshold = 0.02:0.02:0.3; % for reduction

doSave = true;
doPlot = true;
doPrint = true;

% Run through events
for iEv = 1%31%2:numel(iFF_with_DF)
  iFF = iFF_with_DF(iEv);
  str_event = sprintf('Ev_iFF=%04.0f_nMovMean=%g_mask=%.0f_%.0f',iFF,nMovMean,mask_nMovMean(1),mask_nMovMean(1));
  
  % Define time
  t0 = EpochTT(db_table.time(iFF));
  T =  db_table.t_duration_s(iFF);  
  t_ff_start = t0;
  t_ff_stop = t_ff_start + T;
  tint = [t_ff_start t_ff_stop] + 5*[-1 1];
  %tint = tint + [-90 0];
  tDF = int64(db_table.t_df(iFF));
  tDF_EpochTT = EpochTT(tDF);
  
  disp('------------------------------------------------------------')
  if isnan(db_table.t_df(iFF)), strTDF = 'N/A'; else strTDF = EpochTT(int64(db_table.t_df(iFF))); end
  disp("iFF = " + iFF + ...
    ", t1 = " + tint(1).utc('yyyy-mm-ddTHH:MM:SS') + ...
    ", t2 = " + tint(2).utc('yyyy-mm-ddTHH:MM:SS') + ...
    ", tDF = " + strTDF.utc('yyyy-mm-ddTHH:MM:SS'))
  
  
  disp('>> Loading data.')
  try
    c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
    if isempty(scPot), disp('Empty scPot. Skipping.'); continue; end
    c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
    c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
    iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;
    %iPDist_counts.data(iPDist_counts.data == NaN) = 0;
    %iPDist_counts = iPDist_counts.nan2zero;
  catch ME
     warning('Error loading data for ic=%g. Continuing to next event.',ic)
     %continue
  end

  disp('>> Loading B and computing curvature.')
  try
    %%
    c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
    c_eval('gseR? = mms.get_data(''R_gse'',tint+60*[-1 1],?);',1:4)
    c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',1:4)
    c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',1:4)
    c_eval('R? = mms_dsl2gse(gseR?.resample(iPDist),defatt?,-1);',1:4)
    c_eval('B? = dmpaB?.resample(iPDist);',1:4)
    [curvB,avB]=c_4_grad('R?','B?','curvature');
  catch ME
    warning([ME.message ' >> Proceeding without B and curvature.'])
    curvB = TSeries([]);
    avB = TSeries([]);
  end

  %% Prep data
  disp('>> Prepping PDist.')
  PD_use = iPDist.elim([50 Inf]);

  PDc_use = iPDist_counts.nan2zero;

  mat = iPDist_counts.find_low_counts('counts',mask_Counts,'nMovMean',mask_nMovMean,'output','mat');
  
  PDmask = PD_use.mask('energy','mat',mat);  
  PDmm = PDmask.movmean(nMovMean); % see how this affects the results
  PD = PDmm(1:nMovMean:PD_use.length); % This step is done just to reduce the computations for now, need to check if it is needded for noise
  
  PDcmask = PDc_use.mask('energy','mat',mat);  
  PDcmm = PDcmask.movmean(nMovMean);
  PDc = PDcmm(1:nMovMean:PDc_use.length);

  disp('>> Generating microparticles.')
  MP = PD.macroparticles('ntot',nMP,'skipzero',1, ...
    'scpot',scPot.resample(PD),'v_init','v2','MParts', ...
    'counts',PDc);
  MP.userdata.parent_pdist.maskCounts = mask_Counts;
  MP.userdata.parent_pdist.maskMovMean = mask_nMovMean;
  MP.userdata.parent_pdist.movMean = nMovMean;

  disp('>> Reducing distributions.')
  pdist_derived.fx = PD.reduce('1D',[1 0 0]);
  pdist_derived.fy = PD.reduce('1D',[0 1 0]);
  pdist_derived.fz = PD.reduce('1D',[0 0 1]);
  pdist_derived.omni = PD.omni;

  %% Run GMM
  disp('>> Running fit.')
  options = statset('MaxIter',150,'TolFun',1e-5);
  GM = GMDist.fit(MP,K,'start','randSample', ...
    'Options',options,'SharedCovariance',false, ...
    'RegularizationValue',0.1,'Replicates',5);

  %%% Cache merge error
  %fprintf(sprintf('>> Caching merge error, %s.\n',reduction.metric))  
  %GM = GM.cacheMergeErr('metric',reduction.metric);        % the slow step, once  

  %% Save ancillary data to GM object
  fprintf('>> Cache merge error, %s.\n',reduction.metric)
  GM = GM.cacheMergeErr('metric',reduction.metric);        % the slow step, once  

  fprintf('>> Cache B, curvB, n.\n')  
  GM = GM.cacheTS('B',avB); 
  GM = GM.cacheTS('curvB',curvB); 
  GM = GM.cacheTS('n',MP.n); 
  GM = GM.cacheTS('nMP',irf.ts_scalar(GM.time,MP.N));
  
  GM.userdata.pdist_derived = pdist_derived;
  GM.userdata.db_tFF_start = t_ff_start;
  GM.userdata.db_tFF_stop = t_ff_stop;
  GM.userdata.db_tDF = tDF_EpochTT;
  GM.userdata.iSC = ic; 

  GM_save = GM.removeField('data'); % Remove MParts for space

  %% Save GM
  if doSave 
    save([pathSave 'GMDist_' str_event '.mat'],'GM_save'); 
    fprintf('Saved %s.\n',[pathSave 'GMDist_' str_event '.mat'])
  end 

  if doPlot
    %% Figure, OMNIDEFlux, f(vz), BIC, K, M(threshold)
    fig = figure(101);
    fig.Position = [716     1   712   841];
    str_print = sprintf('%s_M_dep_on_threshold',str_event);
  
    nCompPreMerge = GM.nComp;
    nCompMax = max(nCompPreMerge(:));
    nCompMin = min(nCompPreMerge(:));
  
    h = irf_plot(3+nCompMax);
    
    hca = irf_panel('omni fpi');
    irf_spectrogram(hca,PD.deflux.omni.specrec)
  
    hca = irf_panel('BIC'); % it is easier to see counter-streaming populations in this one
    irf_spectrogram(hca,GM.tsBIC('normalization','min').specrec('p_label','BIC-min(BIC)'),'lin')
    
    hca = irf_panel('fz fpi'); % it is easier to see counter-streaming populations in this one
    irf_spectrogram(hca,pdist_derived.fz.specrec)
    %hca.CLim = prctile(hca.Children.CData(:),[1 99]);  
  
    for K_ = nCompPreMerge(1,:)
      G = GM.slice(':',K_);
      hca = irf_panel(sprintf('fz K = %g',K_));
      specrec = G.tsReduced(3).specrec;
      specrec.p_label = {'PSD (s/m^4)';sprintf('K = %g',K_)};
      [~,hcb] = irf_spectrogram(hca,specrec,'log');
      %irf_legend(hca,sprintf('K = %g',K),[0.02 0.1])       
    end
  
    hlinks = linkprop(h(3:end),{'CLim','YLim'});
  
    colormap([flipdim(irf_colormap('Spectral'),1)])
    c_eval('h(?).Layer = ''top'';',1:numel(h))
    irf_plot_axis_align
    irf_zoom(h,'x',PD.time)
    irf_timeaxis(h)
    h(end).XTickLabelRotation = 0;
  
    if doPrint, cn.print(str_print); end 
  end
end

irf.log('on')