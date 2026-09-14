function out = gmm_compare_vdfs_envelope(g1,g2,varargin)
% Idea is to compare the overall end gaussianity, or difference between
% original and merged gmm, to see how many groups are needed for a given
% end threshold on maxwellianity.

% Defaults
method_points = 'montecarlo';
method_comparison = 'rms';
N = 100000;
vvec = -2500:50:2500;
have_options = 0;

% Check inputs
nargs = numel(varargin);
if nargs > 0, have_options = 1; args = varargin(:); end

while have_options
    switch lower(args{1})
      case {'method_points'}
        method_points = args{2};               
        args = args(3:end);
      case {'method_comparison'}
        method_comparison = args{2};               
        args = args(3:end);
      case {'nmc'}
        N = args{2};
        args = args(3:end);
      case {'grid'}
        vvec = args{2};
        args = args(3:end);
      otherwise
        error('Unknown option: %s', args{1});
    end    
  if isempty(args), break, end
end

% Check dimensionality
if iscell(g1)
  [nt, nK] = size(g1);
elseif isa(g1,'gmdistribution')
  nt = 1;
  nK = 1;
end
if iscell(g2)
  nGmax = size(g2,3);
else
  nGmax = 1;
end

f_out = nan(nt,nK,nGmax); % need to change this to cell
for it = 1:nt % Loop over time
  for iK = 1:nK
    if iscell(g1{it,iK})
      nG = numel(g1{it,iK});
    elseif isa(g1{it,iK},'gmdistribution')
      nG = g1{it,iK}.NumComponents;
    end
    for iG = 1:nG
      if iscell(g1{it,iK})
        g1_tmp = g1{it,iK}{iG};
      else
        g1_tmp = g1{it,iK};
      end
      g2_tmp = g2{it,iK,iG};
      if isempty(g1_tmp); disp(sprintf('it = %g, iK = %g, iG = %g: empty gm1',it,iK,iG)); continue; end
      if isempty(g2_tmp); disp(sprintf('it = %g, iK = %g, iG = %g: empty gm2 (merged)',it,iK,iG)); continue; end
      f_tmp = gmm_compare_vdfs(g1_tmp,g2_tmp,varargin{:});
      f_out(it,iK,iG) = f_tmp;
    end
  end
end
out = f_out;



