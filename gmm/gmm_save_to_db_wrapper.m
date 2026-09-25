file_db = strsplit(printpath,'/'); file_db = strrep([strjoin(file_db(1:5),filesep) filesep 'Research/GMM/GMM_BFFs/db_gmm.mat'],'''','"');

% clear DB
DB.metadata.reduction = [];
DB.metadata.list_events = [];
for K = 1:10
  [partitions,unique_groups,map_part2group] = allPartitions(K);
  DB.metadata.reduction(K).numPopulations = K;
  for iGroup = 1:numel(unique_groups)
    DB.metadata.reduction(K).groups(iGroup).id = iGroup;
    DB.metadata.reduction(K).groups(iGroup).members = unique_groups(iGroup);
  end
  for iPartition = 1:numel(partitions)
    DB.metadata.reduction(K).partitions(iPartition).id = iPartition;
    DB.metadata.reduction(K).partitions(iPartition).members = map_part2group(iPartition);
  end
  %DB.metadata_reduction(K).partitions = partitions;
  %DB.metadata_reduction(K).map_partitions_to_groups = map_part2group;
end
save(file_db,'DB')
%%
DB = load(file_db); DB = DB.DB;
id_event = iDF;
id_run = 'run_00'; 

nMacroParticles = cat(1,allMP(:).N);
for iK = 1:size(gm,2)
  id_run = sprintf('run_K=%02.0f',vecK(iK)); 
  % Save metadata into a structure
  event_metadata.time_start = PD(1).time;
  event_metadata.time_stop = PD(PD.length).time;
  event_metadata.nMacroParticles = nMacroParticles; % (nVDF x 1)
  %run_metadata = [];
  
  % Add run to database
  DB = gmm_db_add_run(DB, gm(:,iK), event_metadata, run_metadata);
  %DB = gmm_db_add_run(DB, id_event, id_run, 'gmm',gm, nMacroParticles);
end


save(file_db,'DB')

%% Diagnostic figure
figure(33)
h = irf_plot(4);
id_event = 9;

BIC = cat(2,DB.events(id_event).runs(:).BIC);
Ks = cat(2,DB.events(id_event).runs(:).NumComponents);

nt = size(BIC,1);
BICknee_K = zeros(nt,1);
BICknee_BIC = zeros(nt,1);

BICknee_K_d2 = zeros(nt,1);
BICknee_BIC_d2 = zeros(nt,1);

BICrelim_K = zeros(nt,1);
BICrelim_BIC = zeros(nt,1);

BICreldr_K = zeros(nt,1);
BICreldr_BIC = zeros(nt,1);

BICreldr_K_last = zeros(nt,1);
BICreldr_BIC_last = zeros(nt,1);

relativeImprovement = -diff(BIC,1,2) ./ abs(BIC(:,1:end-1));
relativeImprovementThreshold = 0.001;

dBIC = -diff(BIC,1,2);
relativeDrop = dBIC(:,2:end) ./ dBIC(:,1:end-1);
relativeDropThreshold = 0.4;

for it = 1:nt;
  % Option 1
  [res_x, idx_of_result] = knee_pt(Ks,BIC(it,:));
  BICknee_K(it) = Ks(idx_of_result);
  BICknee_BIC(it) = res_x;

  % Option 2
  d2BIC = diff(BIC(it,:), 2);
  idx = find(d2BIC(1:end-1).*d2BIC(2:end) < 0, 1); 
  if not(isempty(idx))
    BICknee_K_d2(it) = Ks(idx+1);
    BICknee_BIC_d2(it) = BIC(it,idx+1);
  else
    BICknee_K_d2(it) = NaN;
    BICknee_BIC_d2(it) = NaN;
  end

  % Option 3
  [idx] = find(relativeImprovement(it,:)<relativeImprovementThreshold,1,'first');
  if not(isempty(idx))
    BICrelim_K(it) = Ks(idx+0);
    BICrelim_BIC(it) = BIC(it,idx+0);
  else
    BICrelim_K(it) = NaN;
    BICrelim_BIC(it) = NaN;
  end


  % Option 4
  [idx] = find(relativeDrop(it,:)<relativeDropThreshold,1,'first');
  if isempty(idx), idx = 1; end
  BICreldr_K(it) = Ks(idx+1);
  BICreldr_BIC(it) = BIC(it,idx+1);

  % Option 5, also relative drop, but find the last one
  relativeDrop_pos = relativeDrop(it,:);
  relativeDrop_pos(relativeDrop_pos<0) = NaN;
  [idx] = find(relativeDrop_pos<relativeDropThreshold,1,'last');
  if isempty(idx), idx = 1; end
  BICreldr_K_last(it) = Ks(idx+1);
  BICreldr_BIC_last(it) = BIC(it,idx+1);

  if 0
    %%
    %it = 40;
    figure(32)
    h = setup_subplots(3,1); isub = 1;

    hca = h(isub); isub = isub + 1;
    plot(hca,Ks,BIC(it,:),'-*')
    hca.YLabel.String = 'BIC';

    hca = h(isub); isub = isub + 1;
    plot(hca,Ks(1:end-1)+0.5*diff(Ks),dBIC(it,:),'-*')
    hca.YLabel.String = '-\Delta BIC';
    
    hca = h(isub); isub = isub + 1;
    plot(hca,Ks(2:end-1),relativeDrop(it,:),'-*',Ks,Ks*0+relativeDropThreshold)
    hca.YLabel.String = '\Delta BIC(K+1)/\Delta BIC(K)';
    %hca.YLim = [0 2];
    %hca.YTick = 0:0.2:2;
    hca.YGrid = 'on';
    

    hlinks = linkprop(h,{'XLim'});


  end
end
tsBICkneeK = irf.ts_scalar(PD.time,BICknee_K);
tsBICkneeK_d2 = irf.ts_scalar(PD.time,BICknee_K_d2);
tsBICrelimprov = irf.ts_scalar(PD.time,BICrelim_K);
tsBICreldrop = irf.ts_scalar(PD.time,BICreldr_K);

tsBICreldrop_last = irf.ts_scalar(PD.time,BICreldr_K_last);



if 1 % B
  hca = irf_panel('B');
  irf_plot(hca,gseB)
  set(hca,'ColorOrder',mms_colors('xyza'))
  hca.YLabel.String = 'B (nT)';
end
if 1 % dEFlux ion PD2
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
end

if 1 % fred by mms
  hca = irf_panel('fred z mms');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_spectrogram(hca,vdf_fz.specrec,'log')   
  hca.YLabel.String = 'v_z (km/s)';  
  irf_legend(hca,'MMS',[0.98 0.98],'k')
end
if 1 % BIC(K)
  hca = irf_panel('BIC');
  BIC = cat(2,DB.events(id_event).runs(:).BIC);
  BIC = real(BIC);
  maxBIC = repmat(BIC(:,1),[1 size(BIC,2)]);
  specrec.p = BIC./maxBIC;
  specrec.f = cat(2,DB.events(id_event).runs(:).NumComponents);  
  specrec.t = PD.time.epochUnix;
  irf_spectrogram(hca,specrec,'lin')  

  hold(hca,'on')
  %irf_plot(hca,tsBICkneeK,'k')
  %irf_plot(hca,tsBICkneeK_d2,'b')
  %irf_plot(hca,tsBICrelimprov,'r')
  irf_plot(hca,tsBICreldrop,'k')
  %irf_plot(hca,tsBICreldrop_last,'b')
  hold(hca,'off')

  hca.YLabel.String = 'K';
  hcb = colorbar(hca);
  hcb.YLabel.String = 'BIC';
  irf_legend(hca,sprintf('rel. improv. < %.2f',relativeDropThreshold),[0.98 0.1],'fontsize',12,'color','w')
end

c_eval('h(?).Layer = ''top'';',1:numel(h))
colormap([flipdim(irf_colormap('Spectral'),1)])
irf_plot_axis_align(h)
irf_zoom(h,'x',PD.time([1 PD.length]))
h(end).XTickLabelRotation = 0;
%irf_timeaxis(h)

%
if 0
  %%
  cmap = [flipdim(irf_colormap('Spectral'),1)];
  cmap_new = zeros(PD.length,3);
  cmap_new(:,1) = interp1(1:size(cmap,1),cmap(:,1),linspace(1,size(cmap,1),PD.length));
  cmap_new(:,2) = interp1(1:size(cmap,1),cmap(:,2),linspace(1,size(cmap,1),PD.length));
  cmap_new(:,3) = interp1(1:size(cmap,1),cmap(:,3),linspace(1,size(cmap,1),PD.length));
  

  figure(31); 
  h = setup_subplots(3,2); isub = 1;
  
  BIC = specrec.p;
  maxBIC = repmat(BIC(:,1),[1 size(BIC,2)]);
  
  hca = h(isub); isub = isub + 1;
  plot(hca,BIC')
  set(hca,'ColorOrder',cmap_new)
  hold(hca,'on')
  plot(hca,BICknee_K,BICknee_BIC,'*')
  hold(hca,'off')
  set(hca,'ColorOrder',cmap_new)
  hca.XLabel.String = 'K';
  hca.YLabel.String = 'BIC';

  hca = h(isub); isub = isub + 1;
  plot(hca,BIC')
  set(hca,'ColorOrder',cmap_new)
  hold(hca,'on')
  plot(hca,BICrelim_K,BICrelim_BIC,'*')
  hold(hca,'off')
  set(hca,'ColorOrder',cmap_new)
  hca.XLabel.String = 'K';
  hca.YLabel.String = 'BIC';

  hca = h(isub); isub = isub + 1;
  plot(hca,BIC'./maxBIC')
  hca.XLabel.String = 'K';
  hca.YLabel.String = 'BIC-BIC(K=Kmin)';
  set(hca,'ColorOrder',cmap_new)

  hca = h(isub); isub = isub + 1;
  semilogy(hca,relativeImprovement')
  hca.XLabel.String = 'K';
  hca.YLabel.String = '\Delta_K BIC';
  hca.Title.String = 'relativeImprovement';
  set(hca,'ColorOrder',cmap_new)

  hca = h(isub); isub = isub + 1;
  semilogy(hca,relativeDrop')
  hca.XLabel.String = 'K';
  hca.YLabel.String = '\Delta_K BIC';
  hca.Title.String = 'relativeDrop';
  set(hca,'ColorOrder',cmap_new)

end
% How to uniquely identify events? 
%   Start time and stop time in ttns? Based on first and last VDF?
%   From there it can be connected to Louis's database.
% 
% Structure
% DB
% ├── metadata.reduction
% │   └── K
% │       ├── groups
% │       │   ├── id
% │       │   └── members
% │       ├── partitions
% │           ├── id
% │           ├── members % sufficient with group_ids here
% │           └── % map_partitions_to_groups
% └── events
%     ├── event_001
%     │   ├── metadata
%     │   │   ├── event_id
%     │   │   ├── time_start
%     │   │   ├── time_stop
%     │   │   ├── 
%     │   │   ├── 
%     │   │   ├── 
%     │   │   ├── 
%     │   │   └── 
%     │   ├── data
%     │   └── runs
%     │       ├── run_001
%     │       │   ├── settings
%     │       │   ├── GMM
%     │       │   ├── reduction
%     │       │   └── cold
%     │       └── run_002
%     └── event_002