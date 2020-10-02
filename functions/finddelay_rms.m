function varargout = finddelay_rms(s1,s2,varargin)
% FINDDELAY_RMS Find the lag corresponding to the minimized root mean 
%   square between two signals.
%
% [C,LAG] = finddelay_rms(s1,s2);
% [C,LAG] = finddelay_rms(s1,s2,maxlags);
%
% Example:
%   t1 = 0;
%   t2 = 1;
%   T = 0.9;
%   nt = 300;
%   t = linspace(t1,t2,nt);
%   dt = t(2)-t(1);
%   f = 1/T; 
%   tlag = 0.2;
%   s1 = cos(2*pi*t*f)+0.1*randn(1,nt);
%   s2 = cos(2*pi*(t-tlag)*f)+0.1*randn(1,nt);
%   
%   maxlagt = 0.5;
%   maxlags = fix(maxlagt/dt);
%   [C,LAG] = finddelay_rms(s1,s2,maxlags);
%   tlag_corr = LAG(find(abs(C)==min(abs(C))))*dt;
%   hca = subplot(2,1,1);
%   plot(hca,t,s1,t,s2,t+tlag_corr,s2)
%   hca.XLabel.String = 't';
%   hca = subplot(2,1,2);
%   plot(hca,LAG*dt,C)
%   hca.XLabel.String = 'dt';
%   hca.YLabel.String = 'RMS';

ns = numel(s1);

if numel(varargin) > 0
  maxlag = varargin{1};
else
  maxlag = ns-1;
end

rms = zeros(2*maxlag+1,1);

istart1 = [ones(1,ns) 2:ns];
istop1 = [1:ns ns*ones(1,ns-1)];
istart2 = ns+1-istop1;
istop2 = ns+1-istart1;
lag = istart1-istart2;
keep = find(abs(lag)<=maxlag);
%plot((lag(keep)),'*')
istart1 = istart1(keep);
istart2 = istart2(keep);
istop1 = istop1(keep);
istop2 = istop2(keep);
lag = istart1-istart2;

for ii = 1:numel(istop1)
  %ii
  s1_tmp = s1(istart1(ii):istop1(ii));
  s2_tmp = s2(istart2(ii):istop2(ii));
  rms_tmp = sqrt(sum((s1_tmp-s2_tmp).^2))/numel(s1_tmp);
  rms(ii) = rms_tmp;
  lags(ii) = lag(ii);%istart1(ii)-istart2(ii);
end
varargout{1} = rms;
varargout{2} = lags;
  
  