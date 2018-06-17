function [t,y] = get_time(varargin)
% GET_TIME Collects data from plot by mouse click.
%   [t,data] = GET_TIME - collects until return key is pressed
%   [t,data] = GET_TIME(N) - collects N points
%   [t,data] = GET_TIME(N,'zoom') - zooms in after each pair of times is
%                   chosen, saves all in Nx2 array

[ax,args,nargs] = axescheck(varargin{:});
if ~isempty(ax); axis(ax); end
userdata = get(gcf,'userdata');

% Default values
doZoom = 0;
doEpochTT = 0;
N = [];

% Check for input
if isnumeric(args{1})
  N = args{1}; args = args(2:end); % number of times to collect
end 

have_options = nargs > 1;
while have_options
  l = 0;
  switch(lower(args{1}))
    case {'zoom'} % time
      l = 1;
      doZoom = 1;
    case 'epochtt' % number of Monte Carlo iterations
      l = 1;
      doEpochTT = 1;
  end
  args = args((l+1):end);
  if isempty(args), break, end
end

if N > 0
  if doZoom
    N=varargin{1};
    x=zeros(N,2);
    y=zeros(N,2);
    for ii=1:N
      [xprel,yprel] = ginput(2);
      x(ii,1)=xprel(1); x(ii,2)=xprel(2);
      y(ii,1)=yprel(1); y(ii,2)=yprel(2);
      irf_zoom(userdata.subplot_handles,'x',userdata.t_start_epoch + x(ii,:));
    end    
  else
    [x,y] = ginput(N);
  end
else % zoom not implemented for this one
  [x,y] = ginput;
end
% switch nargin
%   case 0
%     [x,y] = ginput;
%   case 1 % N given
%     [x,y] = ginput(varargin{1});
%   case 2 % N and 'zoom' given
%     N=varargin{1};
%     x=zeros(N,2);
%     y=zeros(N,2);
%     for ii=1:N
%       [xprel,yprel] = ginput(2);
%       x(ii,1)=xprel(1); x(ii,2)=xprel(2);
%       y(ii,1)=yprel(1); y(ii,2)=yprel(2);
%       irf_zoom(userdata.subplot_handles,'x',userdata.t_start_epoch + x(ii,:));
%     end
% end

t = userdata.t_start_epoch + x;

if doEpochTT
  t = irf_time(t,'epoch>epochtt');
end
  