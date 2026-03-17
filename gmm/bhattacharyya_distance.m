function out = bhattacharyya_distance(gm, varargin)
%BHATTACHARYYA_DISTANCE Pairwise 3D Bhattacharyya distance between components.
%
% out = bhattacharyya_distance(gm, 'Name', value, ...)
%
% Inputs
%   gm        : either
%               - gmdistribution object, or
%               - cell array of gmdistribution objects (nt x nK) or (nt x 1)
%
% Options (name-value)
%   'sort'        : true/false (default false). If true, sort components by
%                   temperature proxy using gmm_sort() before computing DB.
%   'reg'         : diagonal regularization added if covariances are not SPD
%                   (default 1e-10, will increase if needed).
%   'maxRegTries' : number of times to increase reg by x10 (default 6).
%
% Output (struct)
%   out.DB{it,iK}      : KxK Bhattacharyya distance matrix
%   out.isort{it,iK}   : component sort order used (sorted->unsorted)
%   out.sort           : sort flag
%   out.minDB(it,iK)   : minimum off-diagonal DB
%   out.meanDB(it,iK)  : mean off-diagonal DB
%

% Defaults
doSort = false;
reg0 = 1e-10;
maxRegTries = 6;

if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'sort'
        doSort = varargin{iarg+1};
      case 'reg'
        reg0 = varargin{iarg+1};
      case 'maxregtries'
        maxRegTries = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end

% Normalize input to cell array gm_cell (nt x nK)
if isa(gm, 'gmdistribution')
  gm_cell = {gm};
elseif iscell(gm)
  gm_cell = gm;
else
  error('gm must be a gmdistribution or a cell array of gmdistribution objects.');
end

[nt, nK] = size(gm_cell);
out = struct();
out.sort = doSort;
out.DB = cell(nt, nK);
out.isort = cell(nt, nK);
out.minDB = nan(nt, nK);
out.meanDB = nan(nt, nK);

for it = 1:nt
  for iK = 1:nK
    g = gm_cell{it, iK};
    if isempty(g)
      continue
    end
    if ~isa(g, 'gmdistribution')
      error('gm{%d,%d} is not a gmdistribution.', it, iK);
    end
    K = g.NumComponents;

    if doSort
      isort = gmm_sort(g);
    else
      isort = 1:K;
    end
    out.isort{it,iK} = isort;

    mu = g.mu(isort,:);         % Kx3
    Sigma = g.Sigma(:,:,isort); % 3x3xK

    DB = nan(K,K);
    DB(1:K+1:end) = 0;
    for i = 1:K
      for j = i+1:K
        dmu = (mu(i,:)-mu(j,:))';
        Si = Sigma(:,:,i);
        Sj = Sigma(:,:,j);
        S = 0.5*(Si + Sj);

        reg = 0;
        for tries = 0:maxRegTries
          reg = (tries==0) * 0 + (tries>0) * (reg0 * 10^(tries-1));
          [RS,pS] = chol(S + reg*eye(3));
          [R1,p1] = chol(Si + reg*eye(3));
          [R2,p2] = chol(Sj + reg*eye(3));
          if pS==0 && p1==0 && p2==0
            break
          end
        end
        if pS~=0 || p1~=0 || p2~=0
          error('Failed to make covariances SPD at it=%d iK=%d pair (%d,%d).', it, iK, i, j);
        end

        logdetS = 2*sum(log(diag(RS)));
        logdetSi = 2*sum(log(diag(R1)));
        logdetSj = 2*sum(log(diag(R2)));

        term1 = 0.125 * (dmu' * ((S + reg*eye(3))\dmu));
        term2 = 0.5 * (logdetS - 0.5*(logdetSi + logdetSj));
        db = term1 + term2;

        DB(i,j) = db;
        DB(j,i) = db;
      end
    end

    out.DB{it,iK} = DB;

    tmp = DB;
    tmp(1:K+1:end) = NaN;
    out.minDB(it,iK) = min(tmp,[],'all','omitnan');
    out.meanDB(it,iK) = mean(tmp(~isnan(tmp)),'omitnan');
  end
end

% If input was a single gmdistribution, collapse outputs for convenience
if isa(gm, 'gmdistribution')
  out.DB = out.DB{1,1};
  out.isort = out.isort{1,1};
  out.minDB = out.minDB(1,1);
  out.meanDB = out.meanDB(1,1);
end

end
