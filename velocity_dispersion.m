function varargout = velocity_dispersion(varargin)
% VELOCITY_DISPERSION Obtain distance to common acceleration point through
% velocity dispersion signatures.

% Defaults
getTimes = 0;

% Collect input
nargs = numel(varargin);      
have_options = 0;
if nargs > 0, have_options = 1; args = varargin(:); end
      
while have_options
  l = 0;
  switch(lower(args{1}))     
    case 'click'
      getTimes = 1;
      n_points = args{2};
      l = 2;    
    case 'tv'
      getTimes = 0;
      tt = args{2};
      vv = args{3};
      n_points = tt.length;
      l = 3;
  otherwise
    l = 1;
    irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
  end
  args = args(l+1:end);    
  if isempty(args), break, end    
end
      
if getTimes
  [tt,vv] = get_time(n_points,'EpochTT');
end

tc = tt(1:(n_points-1));
vc = zeros(n_points-1,1);
t0 = zeros(n_points-1,1);
x0 = zeros(n_points-1,1);
dt = zeros(n_points-1,1);

for ii = 1:(n_points-1)      
  dt(ii) = tt(ii+1)-tt(ii);  
  t0(ii) = vv(ii+1)*dt(ii)/(vv(ii+1)-vv(ii));    
  x0(ii) = vv(ii)*t0(ii);
  vc(ii) = (vv(ii+1)+vv(ii))/2;
end
tc = tc + 0.5*dt;
struct.times = tt;
struct.speeds = vv;
struct.tcenter = tc;
struct.vcenter = vc;
struct.dt = dt;
struct.t0 = t0;
struct.x0 = x0;
struct.ts_dispersion = irf.ts_scalar(tt,vv);
struct.ts_t0 = irf.ts_scalar(tc,t0);
struct.ts_x0 = irf.ts_scalar(tc,x0);



varargout{1} = struct;

%     tsIonDisp = irf.ts_scalar(tint_dispersion,ion_speeds);
%     hold(hca,'on')
%     hline = irf_plot(hca,tsIonDisp);
%     hline.Color = [0 0 0];
%     hline.Marker = '.';
%     hline.LineWidth = 2;
%     hold(hca,'off')    