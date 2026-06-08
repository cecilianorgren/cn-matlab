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
  nG = size(g2,3);
else
  nG = 1;
end

f_out = zeros(nt,nK,nG); % need to change this to cell
for it = 1:nt % Loop over time
  for iK = 1:nK
    for iG = 1:nG
      g1_tmp = g1{it,iK};
      g2_tmp = g2{it,iK,iG};
      f_tmp = gmm_compare_vdfs(g1_tmp,g2_tmp,varargin{:});
      f_out(it,iK,iG) = f_tmp;
    end
  end
end
out = f_out;



