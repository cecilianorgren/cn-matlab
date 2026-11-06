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
%                  {'bicmin','bicchord'}). Names: 'bicmin', 'bicchord', 'cold', 'l1'.
%                  A handle @(GM) must return a struct of arrays with GM.nT rows,
%                  one row per VDF, so that new statistics need no other changes.
%     'intervals'  struct of function handles @(GM) -> logical mask [nT x 1] or
%                  time indices (default: all, pre, post DF from GM.userdata)
%     'coldArgs'   arguments to GMDist.cold (default {'maxColdEV',500,'minVabs',50})
%     'thrL1'      thresholds of the 'l1' statistic (default 0.05:0.05:0.3)

p = inputParser;
addParameter(p,'stats',{'bicmin','bicchord'});
addParameter(p,'intervals',defaultIntervals());
addParameter(p,'coldArgs',{'maxColdEV',500,'minVabs',50});
addParameter(p,'thrL1',0.05:0.05:0.3);
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

out = struct;
out.nEvents = numel(listFiles);
out.nVDFs = struct;
for j = 1:numel(intervalNames), out.nVDFs.(intervalNames{j}) = 0; end
for k = 1:nStat
  for j = 1:numel(intervalNames), out.(statField{k}).(intervalNames{j}) = []; end
end

fprintf('iFile = ')
for iFile = 1:numel(listFiles)
  fprintf('%g,',iFile)
  fileName = [listFiles(iFile).folder filesep listFiles(iFile).name];
  tmp = load(fileName,'GM_save');
  GM = tmp.GM_save;

  % Time masks of the intervals, once per event
  masks = struct;
  for j = 1:numel(intervalNames)
    masks.(intervalNames{j}) = toMask(intervals.(intervalNames{j})(GM),GM.nT);
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
end
fprintf('\n')
end

% ---------------------------------------------------------------- intervals
function intervals = defaultIntervals()
  intervals.all  = @(GM) true(GM.nT,1);
  intervals.pre  = @(GM) GM.time.tlim([GM.userdata.db_tFF_start GM.userdata.db_tDF]);
  intervals.post = @(GM) GM.time.tlim([GM.userdata.db_tDF GM.userdata.db_tFF_stop]);
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
    case 'bicmin',   f = @(GM) statBIC(GM,'min');
    case 'bicchord', f = @(GM) statBIC(GM,'chord');
    case 'cold',     f = @(GM) statCold(GM,opt.coldArgs);
    case 'l1',       f = @(GM) statL1(GM,opt.thrL1);
    otherwise
      error('gmmCrossEventStatistics:stat','Unknown statistic ''%s''.',name)
  end
end

function S = statBIC(GM,method)
  % Number of components chosen by BIC (minimum or chord knee), and the
  % temperature and proportion of the coldest component at that choice.
  G = GM.selectNComp(GM.kneeK('method',method));   % K -> column, not kneeK itself
  pop = G.sortBy('T').populations();               % components coldest first
  S.M = G.nComp;                                   % [nT x 1]
  S.minT = pop.Tiso(:,1);                          % eV
  S.wMinT = pop.prop(:,1);
end

function S = statCold(GM,coldArgs)
  % Weight and number of cold components at the BIC-minimum choice.
  G = GM.selectNComp(GM.kneeK('method','min'));
  isCold = G.cold(coldArgs{:});                    % {nT x 1} of logical [1 x K]
  S.fracCold = nan(G.nT,1);
  S.nCold = nan(G.nT,1);
  for it = 1:G.nT
    g = G.gmm{it};
    if isempty(g), continue, end
    S.fracCold(it) = sum(g.ComponentProportion(isCold{it}));
    S.nCold(it) = nnz(isCold{it});
  end
end

function S = statL1(GM,thr)
  % Number of components M kept by the L1 merge criterion for each threshold,
  % for the highest-K fit of each time.
  GM = GM.cacheMergeErr('metric','l1');            % the slow part, once
  ref = size(GM.gmm,2);
  S.M = nan(GM.nT,numel(thr));
  for k = 1:numel(thr)
    Mk = GM.bestPartition(thr(k),'metric','l1');   % [nT x nC]
    S.M(:,k) = Mk(:,ref);
  end
end