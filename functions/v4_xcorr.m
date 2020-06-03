function varargout = v4_xcorr(R1,R2,R3,R4,A1,A2,A3,A4,tint)
% V4_XCORR
% Do automatic cross correlation between fields
% out = v4_xcorr(R1,R2,R3,R4,A1,A2,A3,A4,tint)

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

% Cross-correlation between individual spacecraft tofind time delay
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
    tic;[corr,lags] = xcorr(D1,D2,'coeff',fix(nt/10)); toc
    [cmax,icmax] = max(corr);
    dtcmax = -lags(icmax)*dt;
    %tdelay{icomp}(isc) = finddelay(D1,D2)*dt;        
    tdelay{icomp}(isc) = dtcmax;
    corrall{icomp}(isc) = cmax;
    
    tic;[corr_rms,lags_rms] = finddelay_rms(D1,D2,fix(nt/10)); toc
    [cmax_rms,icmax_rms] = min(corr_rms);
    dtcmax_rms = -lags_rms(icmax_rms)*dt;
    tdelay_rms{icomp}(isc) = dtcmax_rms;
    corrall_rms{icomp}(isc) = cmax_rms;
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
  
  t_rms = tint(1) + tdelay_rms{icomp};
  v_rms{icomp}  = mms.mms4_v(t_rms);
end
varargout{1} = tdelay;
varargout{2} = v;
varargout{3} = corrall;
varargout{4} = tdelay_rms;
varargout{5} = v_rms;
varargout{6} = corrall_rms;
