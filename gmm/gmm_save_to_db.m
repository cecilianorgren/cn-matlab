function db = gmm_save_to_db(db,id_event,desc_run,type,gmm,nMacroParticles)


%if ~isfield(db.events)
%  db.events = [];
%end

% Check db for the event
if isfield(db.metadata,'events') && ~intersect(db.metadata.events,id_event)
  % Add to list of events
  db.metadata.list_events = sort([db.metadata.list_events id_event]);
  % Add empty event to structure
  db.events(id_event) = [];
end

switch type
  case 'gmm'
    % Check if event exists, otherwise initiate it

    % Save event metadata
    db.events(id_event).metadata.id = id_event;

    % Check event for the run
    if ~isfield(db.events(id_event),'list_runs_id') % new event with no runs
      db.events(id_event).list_runs_id = [];
      db.events(id_event).list_runs_desc = {};
      db.events(id_event).runs = [];
      % Add run
      id_run = 1;
      db.events(id_event).list_runs_id(id_run) = id_run;
      db.events(id_event).list_runs_desc{id_run} = desc_run;
    elseif any(cellfun(@(x)strcmp(x,desc_run),db.events(id_event).list_runs_desc)) % exiting event and existing run
      id_run = find(cellfun(@(x)strcmp(x,desc_run),db.events(id_event).list_runs_desc));   
      db.events(id_event).list_runs_id(id_run) = id_run;
      db.events(id_event).list_runs_desc{id_run} = desc_run;        
    else  % existing event but new run
      id_run = numel(db.events(id_event).list_runs_id) + 1;
      db.events(id_event).list_runs_id(id_run) = id_run;
      db.events(id_event).list_runs_desc{id_run} = desc_run;        
    end

    % Save run metadata
    settings.nMacroParticles = nMacroParticles;
    settings.NumComponents = gmm{1}.NumComponents;
    settings.CovarianceType = gmm{1}.CovarianceType;

    nt = size(gmm,1);
    K = settings.NumComponents;
    weight = zeros(nt,K);
    mu = zeros(nt,K,3);
    Sigma = zeros(nt,3,3,K);
    for it = 1:nt
      weight(it,:) = gmm{it}.ComponentProportion;
      mu(it,:,:) = gmm{it}.mu;
      Sigma(it,:,:,:) = gmm{it}.Sigma;
    end
   
    db.events(id_event).runs(id_run).settings = settings;
    db.events(id_event).runs(id_run).weight = weight;
    db.events(id_event).runs(id_run).mu = mu;
    db.events(id_event).runs(id_run).Sigma = Sigma;
    

end