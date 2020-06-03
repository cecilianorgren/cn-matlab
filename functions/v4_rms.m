function varargout = v4_rms(R1,R2,R3,R4,A1,A2,A3,A4,tint,dtmax)
% V4_RMS
% Calculate the tomedelay that give the minimal rms difference
% out = V4_RMS(R1,R2,R3,R4,A1,A2,A3,A4,tint)

ncomp = size(A1.data,2);
R.time = R1.time;
R.gseR1 = R1;
R.gseR2 = R2;
R.gseR3 = R3;
R.gseR4 = R4;
if isempty(tint)
  tint = A1.time([1 end]);
end
A1 = A1.tlim(tint);
dt = A1.time(2)-A1.time(1);
nt = A1.length;

if isempty(dtmax)
  maxlags = fix(nt/10);
else
  maxlags = fix(dtmax/dt);
end

% Calculate rms difference between individual spacecraft to find time delay
for isc = 1:4
  
  c_eval('B = A?.resample(A1).tlim(tint);',isc)
  for icomp = 1:ncomp
    if ncomp == 3
      comp_str = A1.representation{1}{icomp};                        
      D1 = A1.(comp_str).data;
      D2 = B.(comp_str).data;
    else
      D1 = A1.data;
      D2 = B.data;      
    end
    
    if isc == 1 % sc1 is used as reference, so time delay and rms difference are both zero and no need to calculate
      tdelay{icomp}(isc) = 0;
      corrall{icomp}(isc) = 0;
      continue; % jump tp sc2
    end
  
    tic;[corr,lags] = findlag_rms(D1,D2,maxlags); toc
    [cmax,icmax] = min(corr);
    dtcmax= -lags(icmax)*dt;
    tdelay{icomp}(isc) = dtcmax;
    corrall{icomp}(isc) = cmax;
  end
end
% for icomp = 1:3
%   irf_plot({A1.(comp_str),B.(comp_str)},'comp','dt',tdelay{icomp})
%   pause
% end
%tdelay

% Velocity from time delay and positions
%   Inputs:
%   t  = [t1, t2, t3, t4]; in EpochTT
%   v  = TSeries of [vx, vy, vz]
%     ie: v = irf.ts_vec_xyz(t, [vx, vy, vz]); t in EpochTT and v in GSE ref frame.
%   dt = [0, t2-t1, t3-t1, t4-t1];
for icomp = 1:ncomp
  t = tint(1) + tdelay{icomp};
  v{icomp}  = mms.mms4_v(t);  
end
varargout{1} = tdelay;
varargout{2} = v;
varargout{3} = corrall;


end
function varargout = findlag_rms(s1,s2,varargin)
  ns = numel(s1);

  if numel(varargin) > 0
    maxlag = varargin{1};
  else
    maxlag = ns-1;
  end

  rms = zeros(2*maxlag+1,1);
  lags = zeros(2*maxlag+1,1);
  
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
    lags(ii) = istart1(ii)-istart2(ii);    
  end
  plot(lags,rms)
  varargout{1} = rms;
  varargout{2} = lags;
end