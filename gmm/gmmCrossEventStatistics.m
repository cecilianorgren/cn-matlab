function out = gmmCrossEventStatistics(listFiles,varargin)
% GMMCROSSEVENTSTATISTICS Per-VDF statistics of GMDist objects over many events,
%   split into time intervals.
%
%   out = gmmCrossEventStatistics(listFiles[,option,value,...])
%
%   For each event (file with GM_save) every statistic is computed once on the
%   whole event, as one row per VDF, and the rows are then split into the
%   intervals (all, pre, post, ...) and appended to what has been collected.
%
%   Output: out.<statistic>.<interval>.<field> with one row per VDF, plus the
%   bookkeeping fields event (index in listFiles) and iVDF (index in the event).
%   out.nVDFs.<interval> and out.nEvents count what went in.
%
%   Options:
%     'stats'      cell of statistic names or function handles (default
%                  {'bicmin','bicchord'}). Names: 'bicmin', 'bicchord', 'l1'.
%                  Every selection reports the same fields: M, minT, wMinT and, with
%                  coldArgs, nCold and fracCold.
%                  A handle @(GM) must return a struct of arrays with GM.nT rows,
%                  one row per VDF, so that new statistics need no other changes.
%     'intervals'  struct of function handles @(GM) -> logical mask [nT x 1] or
%                  time indices (default: all, pre, post DF from GM.userdata)
%     'coldArgs'   arguments to GMDist.cold (default {'maxColdEV',500,'minVabs',50}),
%                  applied to the selected mixture of each statistic: nCold is the number
%                  of cold components, fracCold their share of the weight; [] skips it
%     'thrL1'      thresholds of the 'l1' statistic (default 0.05:0.05:0.3); the
%                  column k of out.l1.<interval>.M belongs to thrL1(k), also in out.thrL1
%     'l1Args'     extra arguments for the L1 merge errors, e.g. {'nSamples',20000,'seed',1},
%                  to match the stored errors (see GMDist.cacheMergeErr)
%   The 'l1' statistic uses the L1 errors stored in the GMDist (cacheMergeErr);
%   it warns if there are none that match, since they would then be recomputed (slow).
%
%   Selecting data with the TSeries stored in the GMDist (GMDist.cacheTS):
%     'eventFilter'  @(GM) -> logical scalar; events for which it is false are
%                    skipped altogether, e.g. @(GM) mean(GM.getTS('B').data(:,1)) < thrBx
%     'vdfFilter'    @(GM) -> logical mask [nT x 1] (or indices); VDFs outside it
%                    are dropped in all intervals, e.g. @(GM) GM.getTS('B').data(:,1) < thrBx
%   Values to return alongside the statistics ('covariates'): a cell of names of
%   stored TSeries, e.g. {'B','curvB'}, or a struct of handles @(GM) -> [nT x k]
%   (e.g. struct('Bmag',@(GM) vecnorm(GM.getTS('B').data,2,2))). They are
%   resampled to the VDF times (GM.getTS) and returned in out.cov.<interval>.<name>,
%   row for row aligned with out.<statistic>.<interval>, which have the same event
%   and iVDF columns. For example scatter(out.cov.pre.curvB,out.bicmin.pre.minT).
%   out.events(i) tells which events were used (file, used, nVDFs); out.nEventsUsed.

p = inputParser;
addParameter(p,'stats',{'bicmin','bicchord'});
addParameter(p,'intervals',defaultIntervals());
addParameter(p,'coldArgs',{'maxColdEV',500,'minVabs',50});
addParameter(p,'thrL1',0.05:0.05:0.3);
addParameter(p,'l1Args',{});
addParameter(p,'eventFilter',[]);
addParameter(p,'vdfFilter',[]);
addParameter(p,'covariates',{});
parse(p,varargin{:});

statNames = p.Results.stats;
if ~iscell(statNames), statNames = {statNames}; end
intervals = p.Results.intervals;
intervalNames = fieldnames(intervals);

% One function handle per statistic, and the name to store it under
nStat = numel(statNames);
statFun = cell(1,nStat);
statField = cell(1,nStat);
for k = 1:nStat
  if isa(statNames{k},'function_handle')
    statFun{k} = statNames{k};
    statField{k} = sprintf('stat%d',k);
  else
    statFun{k} = statRegistry(statNames{k},p.Results);
    statField{k} = lower(statNames{k});
  end
end

covFun = covariateFunctions(p.Results.covariates);
covNames = fieldnames(covFun);
eventFilter = p.Results.eventFilter;
vdfFilter = p.Results.vdfFilter;

out = struct;
out.nEvents = numel(listFiles);
out.thrL1 = p.Results.thrL1;
out.nVDFs = struct;
for j = 1:numel(intervalNames), out.nVDFs.(intervalNames{j}) = 0; end
out.events = struct('file',{},'used',{},'nVDFs',{});
for j = 1:numel(intervalNames), out.cov.(intervalNames{j}) = []; end
for k = 1:nStat
  for j = 1:numel(intervalNames), out.(statField{k}).(intervalNames{j}) = []; end
end

fprintf('iFile = ')
for iFile = 1:numel(listFiles)
  fprintf('%g,',iFile)
  fileName = [listFiles(iFile).folder filesep listFiles(iFile).name];
  tmp = load(fileName,'GM_save');
  GM = tmp.GM_save;

  % Events to skip altogether
  used = true;
  if ~isempty(eventFilter), used = logical(eventFilter(GM)); end
  out.events(iFile) = struct('file',fileName,'used',used,'nVDFs',GM.nT);
  if ~used, continue, end

  % VDF filter, and the time masks of the intervals, once per event
  keep = true(GM.nT,1);
  if ~isempty(vdfFilter), keep = toMask(vdfFilter(GM),GM.nT); end
  masks = struct;
  for j = 1:numel(intervalNames)
    masks.(intervalNames{j}) = toMask(intervals.(intervalNames{j})(GM),GM.nT) & keep;
    out.nVDFs.(intervalNames{j}) = out.nVDFs.(intervalNames{j}) + nnz(masks.(intervalNames{j}));
  end

  % Each statistic is computed once on the whole event, then split by interval
  for k = 1:nStat
    S = statFun{k}(GM);
    S.event = repmat(iFile,GM.nT,1);
    S.iVDF = (1:GM.nT).';
    for j = 1:numel(intervalNames)
      name = intervalNames{j};
      out.(statField{k}).(name) = appendRows(out.(statField{k}).(name),takeRows(S,masks.(name)));
    end
  end

  % Covariates: the same rows as the statistics
  if ~isempty(covNames)
    C = struct;
    for k = 1:numel(covNames)
      v = covFun.(covNames{k})(GM);
      if size(v,1) ~= GM.nT
        error('gmmCrossEventStatistics:cov','Covariate ''%s'' must have one row per VDF (%d), not %d.',covNames{k},GM.nT,size(v,1))
      end
      C.(covNames{k}) = v;
    end
    C.event = repmat(iFile,GM.nT,1);
    C.iVDF = (1:GM.nT).';
    for j = 1:numel(intervalNames)
      name = intervalNames{j};
      out.cov.(name) = appendRows(out.cov.(name),takeRows(C,masks.(name)));
    end
  end
end
out.nEventsUsed = sum([out.events.used]);
fprintf('\n')
end

% ---------------------------------------------------------------- intervals
function intervals = defaultIntervals()
  intervals.all  = @(GM) true(GM.nT,1);
  intervals.pre  = @(GM) GM.time.tlim([GM.userdata.db_tFF_start GM.userdata.db_tDF]);
  intervals.post = @(GM) GM.time.tlim([GM.userdata.db_tDF GM.userdata.db_tFF_stop]);
  intervals.pre5s  = @(GM) GM.time.tlim(GM.userdata.db_tDF + [-5 -1]);
  intervals.post5s = @(GM) GM.time.tlim(GM.userdata.db_tDF + [1 +5]);
end

function m = toMask(idx,n)
  % Logical mask [n x 1] from a mask or a vector of indices.
  if islogical(idx)
    m = idx(:);
  else
    m = false(n,1);
    m(idx) = true;
  end
end

% ---------------------------------------------------------------- covariates
function f = covariateFunctions(cov)
  % Struct of handles @(GM) -> [nT x k] from names of stored TSeries or a struct.
  f = struct;
  if isempty(cov), return, end
  if isstruct(cov), f = cov; return, end
  if ischar(cov) || isstring(cov), cov = {char(cov)}; end
  for k = 1:numel(cov)
    name = cov{k};
    f.(name) = @(GM) GM.getTS(name).data;           % stored TSeries, resampled to the VDF times
  end
end

% ---------------------------------------------------------------- bookkeeping
function S = takeRows(S,mask)
  f = fieldnames(S);
  for k = 1:numel(f), S.(f{k}) = S.(f{k})(mask,:); end
end

function A = appendRows(A,S)
  if isempty(A), A = S; return, end
  f = fieldnames(S);
  for k = 1:numel(f), A.(f{k}) = [A.(f{k}); S.(f{k})]; end
end

% ---------------------------------------------------------------- statistics
function f = statRegistry(name,opt)
  % Add a statistic here (or pass a handle in 'stats'): it takes one GMDist and
  % returns a struct of arrays with GM.nT rows.
  switch lower(name)
    case 'bicmin',   f = @(GM) statBIC(GM,'min',opt.coldArgs);
    case 'bicchord', f = @(GM) statBIC(GM,'chord',opt.coldArgs);
    case 'l1',       f = @(GM) statL1(GM,opt.thrL1,opt.l1Args,opt.coldArgs);
    otherwise
      error('gmmCrossEventStatistics:stat','Unknown statistic ''%s''.',name)
  end
end

function S = mixtureStats(G,coldArgs)
  % The statistics of a single-column GMDist, one row per time: number of
  % components, temperature and proportion of the coldest component, and (if
  % coldArgs is not empty) the number of cold components and their weight.
  % Every selection (BIC, L1 merge, ...) reports through this, so they give the
  % same fields.
  pop = G.sortBy('T').populations();               % components coldest first
  S.M = G.nComp;                                   % [nT x 1]
  S.minT = pop.Tiso(:,1);                          % eV  
  S.maxT = max(pop.Tiso(:,:),[],2);                % eV
  S.wMinT = pop.prop(:,1);
  S.wMaxT = max(pop.prop(:,:),[],2);
  S.Eth = pop.Eth(:,1);                            % eV
  S.Ebulk = pop.Ebulk(:,1);                        % eV
  S.Etot = pop.Etot(:,1);                          % eV
  if ~isempty(coldArgs)
    isCold = G.cold(coldArgs{:});                  % {nT x 1} of logical [1 x K]
    S.nCold = nan(G.nT,1);
    S.fracCold = nan(G.nT,1);
    for it = 1:G.nT
      g = G.gmm{it};
      if isempty(g), continue, end
      S.nCold(it) = nnz(isCold{it});
      S.fracCold(it) = sum(g.ComponentProportion(isCold{it}));
    end
  end
end

function S = statBIC(GM,method,coldArgs)
  % Components chosen by BIC (minimum or chord knee).
  G = GM.selectNComp(GM.kneeK('method',method));   % K -> column, not kneeK itself
  S = mixtureStats(G,coldArgs);
end

function S = statL1(GM,thr,l1Args,coldArgs)
  % The same statistics as the BIC ones, after merging with the L1 criterion for
  % each threshold (columns of the fields = thresholds), for the highest-K fit of
  % each time, from the L1 errors stored in GM.
  if ~isfield(GM.misc,'partitions') || ~any(strcmp({GM.misc.partitions.metric},'l1'))
    warning('gmmCrossEventStatistics:l1','No stored L1 merge errors in this GMDist: they will be computed (slow).')
  end
  ref = size(GM.gmm,2);
  for k = 1:numel(thr)
    mG = GM.merge('threshold',thr(k),'metric','l1',l1Args{:},'bic',false);   % uses the stored errors
    Sk = mixtureStats(mG.slice(':',ref),coldArgs);
    if k == 1
      S = structfun(@(x) nan(numel(x),numel(thr)),Sk,'UniformOutput',false);
    end
    f = fieldnames(Sk);
    for j = 1:numel(f), S.(f{j})(:,k) = Sk.(f{j}); end
  end
end