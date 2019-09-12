function varargout = get_map_liouville_moments(varargin)

if isempty(varargin)  
  load('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_map_moments.mat');    
elseif isa(varargin{1},'char')
  if strcmp(varargin{1},'load_from_disc')
%     data = load('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_map_moments.mat');
%     fields = fieldnames(data);
%     for ifield = 1:numel(fields)
%       eval(sprintf('%s = data.%s;',fields{ifield},fields{ifield}))
%     end    
    load('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_map_moments.mat');    
  else
    if isa(varargin{1},'char')
      error(sprintf('Input not recognized: %s',varargin{1}))
    else
      disp('Nothing loaded.')
      return;
    end
  end
else
  % function 
  units = irf_units;
  n0 = 0.1;
  T0 = 100; 
  vt0 = sqrt(2*units.e*T0./units.me); % m/s
  v0 = 0;
  n_psi = 189; min_psi = 0.01; max_psi = 100*T0;
  n_vpsi = 129; min_vpsi = 0.01; max_vpsi = 10*vt0*1e-3; % km/s
  doLogspace = 1;
  if doLogspace
    psi = logspace(log10(min_psi),log10(max_psi),n_psi);
    vpsi = logspace(log10(min_vpsi),log10(max_vpsi),n_vpsi);
  else  
    psi = linspace(min_psi,max_psi,n_psi);
    vpsi = linspace(min_vpsi,max_vpsi,n_vpsi);
  end

  [PSI,VPSI] = ndgrid(psi,vpsi);

  n_lb_map = zeros(n_psi,n_vpsi);
  n_sep_map = zeros(n_psi,n_vpsi);
  v_lb_map = zeros(n_psi,n_vpsi);
  v_sep_map = zeros(n_psi,n_vpsi);

  tic;
  for i_psi = 1:n_psi
    fprintf('%g,',i_psi)
    for i_vpsi = 1:n_vpsi
      tmp_psi = psi(i_psi);
      tmp_vpsi = vpsi(i_vpsi);
      [n_lb_tmp,n_sep_tmp,v_lb_tmp,v_sep_tmp] = paper_electron_acceleration.liouville_mapped_nf(n0,T0,v0,tmp_psi,tmp_vpsi);
      n_lb_map(i_psi,i_vpsi) = n_lb_tmp;
      n_sep_map(i_psi,i_vpsi) = n_sep_tmp;
      v_lb_map(i_psi,i_vpsi) = v_lb_tmp;
      v_sep_map(i_psi,i_vpsi) = v_sep_tmp;
    end
  end
  fprintf('\n')
  toc
end
  

% Assign output, regardless of how data was loaded
varargout{1} = PSI/T0;
varargout{2} = VPSI*1e3/vt0;
varargout{3} = n_sep_map*1e-6/n0;
varargout{4} = abs(v_sep_map)/vt0;

