function out = gmm_compare_vdfs(g1,g2,varargin)
% Idea is to compare the overall end gaussianity, or difference between
% original and merged gmm, to see how many groups are needed for a given
% end threshold.

% Defaults
method_points = 'montecarlo';
method_comparison = 'rms';
N = 100000;
vvec = -2500:50:2500;
have_options = 0;

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


switch lower(method_points)
  case 'montecarlo'
    xyz =  gmm_monte_carlo_sampling(g1,N,1:g1.NumComponents);
  case 'grid'
    [X Y Z] = ndgrid(vvec,vvec,vvec);
    xyz = [X(:) Y(:) Z(:)];
end

[f1] = gmm_evaluate_f(g1,xyz,1); % evaluate f at xyz
[f2] = gmm_evaluate_f(g2,xyz,1); % evaluate f at xyz


switch lower(method_comparison)
  case 'rms'
    d = sum(abs(f2-f1)./sum(abs(f1)))/1;
    %d = sum(sqrt((f2-f1).^2)./sqrt((f1).^2))/N;
  case 'epsilon'
    d = [];
end

out = d;



