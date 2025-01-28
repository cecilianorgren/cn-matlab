function varargout = edi_eh_running_average_correlation(A,B,T,dt,maxLag)

B = B.resample(A);

centerTime = B.time(1):dt:B.time(B.length);
nTime = centerTime.length;
vecC = zeros(nTime,1);
vecLag = zeros(nTime,1);

for iT = 1:nTime
  tint = centerTime(iT) + 0.5*T*[-1 1];
  tmpA = A.tlim(tint).data;
  tmpB = B.tlim(tint).data;
  [C,lag] = xcorr(tmpA,tmpB,maxLag,'coeff');
  
  iMaxC = find(C==max(C));
  vecLag(iT) = lag(iMaxC);
  vecC(iT) = C(iMaxC);
  
end

tsC = irf.ts_scalar(centerTime,vecC);
tsC.name = 'Cross correlation';

tsLag = irf.ts_scalar(centerTime,vecLag);
tsLag.name = 'Lag';

varargout{1} = tsC;
if nargout == 2
  varargout{2} = tsLag;
end