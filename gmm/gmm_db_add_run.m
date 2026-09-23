function db = gmm_db_add_run(db,gmm,event_metadata,run_metadata)

% Event unique identifier: time in ns of first and last VDF
t_start = event_metadata.time_start.ttns;
t_stop = event_metadata.time_stop.ttns;


% Check if event already exists
if isfield(db,'events')
  if isfield(db.events,'time_start')
    start_times = cat(1,db.events.time_start);
    stop_times = cat(1,db.events.time_stop);
    [isEvent,id_event] = ismember([t_start t_stop],[start_times stop_times],"rows");
    if ~isEvent % new event      
      nEvents = numel(start_times);
      id_event = nEvents + 1; % incremental id number
    end
  else
    isEvent = false;
    id_event = 1;
  end
else
  isEvent = false;
  id_event = 1;
end 

% Check if run exists within event
if ~isEvent % new event, run will not exist
  id_run = 1;
else % existing event
  if isfield(db.events(id_event),'runs') 
    % How to check if the run is the same...
    % For now, just add runs
    nRuns = numel(db.events(id_event).runs);
    id_run = nRuns + 1;    
  else
    id_run = 1;
  end
end

% 
% % Check event for the run
% if ~isfield(db.events(id_event),'list_runs_id') % new event with no runs
%   db.events(id_event).list_runs_id = [];
%   db.events(id_event).list_runs_desc = {};
%   db.events(id_event).runs = [];
%   % Add run
%   id_run = 1;
%   db.events(id_event).list_runs_id(id_run) = id_run;
%   db.events(id_event).list_runs_desc{id_run} = desc_run;
% elseif any(cellfun(@(x)strcmp(x,desc_run),db.events(id_event).list_runs_desc)) % exiting event and existing run
%   id_run = find(cellfun(@(x)strcmp(x,desc_run),db.events(id_event).list_runs_desc));   
%   db.events(id_event).list_runs_id(id_run) = id_run;
%   db.events(id_event).list_runs_desc{id_run} = desc_run;        
% else  % existing event but new run
%   id_run = numel(db.events(id_event).list_runs_id) + 1;
%   db.events(id_event).list_runs_id(id_run) = id_run;
%   db.events(id_event).list_runs_desc{id_run} = desc_run;        
% end



% Save event metadata
if ~isEvent
  db.events(id_event).time_start = t_start;
  db.events(id_event).time_stop = t_stop;
end


% db.events(id_event).metadata.id = id_event;

% Save run metadata
%settings.nMacroParticles = nMacroParticles;
settings.NumComponents = gmm{1}.NumComponents;
%settings.CovarianceType = gmm{1}.CovarianceType;

nt = size(gmm,1);
K = settings.NumComponents;

% Initialize matrices
weight = zeros(nt,K);
mu = zeros(nt,K,3);
Sigma = zeros(nt,3,3,K);
BIC = zeros(nt,1);
AIC = zeros(nt,1);
NLL = zeros(nt,1);

for it = 1:nt
  weight(it,:) = gmm{it}.ComponentProportion;
  mu(it,:,:) = gmm{it}.mu;
  Sigma(it,:,:,:) = gmm{it}.Sigma;
  BIC(it,1) = gmm{it}.BIC;
  AIC(it,1) = gmm{it}.AIC;
  NLL(it,1) = gmm{it}.NegativeLogLikelihood;
end

db.events(id_event).runs(id_run).NumComponents = K;
db.events(id_event).runs(id_run).weight = weight;
db.events(id_event).runs(id_run).mu = mu;
db.events(id_event).runs(id_run).Sigma = Sigma;
db.events(id_event).runs(id_run).AIC = AIC;
db.events(id_event).runs(id_run).BIC = BIC;
db.events(id_event).runs(id_run).NegativeLogLikelihood = NLL;
db.events(id_event).runs(id_run).gmdistribution = gmm{it};
