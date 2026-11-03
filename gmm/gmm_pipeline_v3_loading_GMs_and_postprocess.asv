%% Load DB table
localuser = 'cecilianorgren';
%file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2022.nc'];


%% Loop events and plot brief overview
units = irf_units;

mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');

doPrint = 1;
doPlot = 0;
vecK = 1:5;

% Turn off irf._log warning messages
irf.log('off')

% Save path
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v1/'];
listFiles = dir([pathSave '*Ev*']);
nFiles = numel(listFiles);

for iFile = 1%:nFiles
  fileName = [listFiles(iFile).folder filesep listFiles(iFile).name];
  GM = load(fileName,'GM_save').GM_save; 

  GBIC = GM.slice(':',GM.kneeK);

  %% Plot, BIC
  h = irf_plot(5);

  hca = irf_panel('omni fpi');
  irf_spectrogram(hca,GM.userdata.pdist_derived.omni.deflux.specrec)

  hca = irf_panel('fz fpi');
  irf_spectrogram(hca,GM.userdata.pdist_derived.fz.specrec)

  hca = irf_panel('BIC/maxBIC)');
  irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec( ...
  'p_label','log BIC/max(BIC)','f_label','K'),'log')
  hold(hca,'on')
  irf_plot(hca,GM.tsBICKnee,'k')
  hold(hca,'off')

  hca = irf_panel('BIC/minBIC');
  irf_spectrogram(hca,GM.tsBIC('normalization','min').specrec( ...
  'p_label','log BIC-min(BIC)','f_label','K'),'log')
  hold(hca,'on')
  irf_plot(hca,GM.tsBICKnee('method','min'),'k')
  irf_plot(hca,GM.tsBICKnee('method','chord'),'r')
  hold(hca,'off')

  hca = irf_panel('fz by BIC');
  irf_spectrogram(hca,GBIC.tsReduced(1).specrec,'log')


  irf_plot_axis_align
  irf_timeaxis(h)
  %% Plot 1

  h = setup_subplots(2,6,'horizontal');
  isub = 1;
  for iK = 1:K
    hca = h(isub); isub = isub + 1;
    G = GM.slice(':',iK).sortBy('T');
    tsN = G.tsN('n',GM.userdata.TSeries.n.data);
    tsT = G.tsT('type','iso');
  
    cellN = cellfun(@(x)x.data,tsN,'UniformOutput',false);
    matN = cat(2,cellN{:});
    cellT = cellfun(@(x)x.data,tsT,'UniformOutput',false);
    matT = cat(2,cellT{:});
    
    scatter(hca,matN,matT,'o')
    hca.XLabel.String = 'n (cm^{-3})';
    hca.YLabel.String = 'T (eV)';
  end  

  for iK = 1:K
    hca = h(isub); isub = isub + 1;
    G = GM.slice(':',iK).sortBy('T');
    tsV = G.tsV;

    cX = cellfun(@(x)x.x.data,tsV,'UniformOutput',false);
    mX = cat(2,cX{:});
    cY = cellfun(@(x)x.z.data,tsV,'UniformOutput',false);
    mY = cat(2,cY{:});
    scatter(hca,mX,mY,'o')
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_z (km/s)';
  end  
  
  c_eval('h(?).Box = ''on'';',1:numel(h))
  c_eval('h(?).FontSize = 12;',1:numel(h))

  %% Plot 2, using .byNComp
  h = setup_subplots(2,6,'vertical');
  isub = 1;

  for iK = 1:K
    G = GM.slice(':',iK).sortBy('T');
    data = G.byNComp('n',GM.userdata.TSeries.n.data);

    hca = h(isub); isub = isub + 1;
    scatter(hca,data.n,data.T,'o')
    %hca.XScale = 'log';
    %hca.YScale = 'log';

    hca = h(isub); isub = isub + 1;
    scatter(hca,data.V(:,:,1),data.V(:,:,3),'o')
  end

  c_eval('h(?).Box = ''on'';',1:numel(h))
  c_eval('h(?).FontSize = 12;',1:numel(h))
end
%%




% Run through events
for iDF = 17:nDF
  iFF = iFFs(iDF);
  str_event = sprintf('Ev_iDF=%03.0f_nMovMean=%g_mask=%.0f_%.0f',iDF,nMovMean,mask_nMovMean(1),mask_nMovMean(1));
  
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
    if isempty(scPot), disp('Empty scPot. Skipping.'); continue; end
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

  MP = PD.macroparticles('ntot',nMP,'skipzero',1, ...
    'scpot',scPot.resample(PD),'v_init','v2','MParts', ...
    'counts',PDc);
  MP.userdata.parent_pdist.maskCounts = mask_Counts;
  MP.userdata.parent_pdist.maskMovMean = mask_nMovMean;
  MP.userdata.parent_pdist.movMean = nMovMean;

  pdist_derived.fx = PD.reduce('1D',[1 0 0]);
  pdist_derived.fy = PD.reduce('1D',[0 1 0]);
  pdist_derived.fz = PD.reduce('1D',[0 0 1]);
  pdist_derived.omni = PD.omni;


  %% Run GMM
  disp('>> Running fit.')
  options = statset('MaxIter',150,'TolFun',1e-5);
  GM = GMDist.fit(MP,K,'start','previous', ...
    'Options',options,'SharedCovariance',false, ...
    'RegularizationValue',0.1);

  %% Reduce mixture
  %disp('>> Reducing mixture.') 
  fprintf(sprintf('>> Caching merge error, %s.\n',reduction.metric))  
  GM = GM.cacheMergeErr('metric',reduction.metric);        % the slow step, once  

  %% Save GM
  GM_save = GM.removeField('data');
  GM_save.userdata.pdist_derived = pdist_derived;
  save([pathSave 'GMDist_' str_event '.mat'],'GM_save')

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
  irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec('p_label','BIC/max(BIC)'),'lin')
  
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

irf.log('on')