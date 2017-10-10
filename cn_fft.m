function [outspecrec,outPxx,outF] = cn_fft(data)
%cn_fft  compute power spectrum for a single time
%
% [t,power,f] = irf_powerfft(data,nfft,sfreq,[overlap])
% [specrec] = irf_powerfft(data,nfft,sfreq,[overlap])
%	SPECREC is a structure:
%		SPECREC.T - time
%		SPECREC.P - spectrum
%		SPECREC.F - frequency
%
% Input data may be either a vector (with first column as time) or a
% TSeries object.
% Please note: Output is still the same, [specrec] struct or [t, power, f]
% with time being of type "double" in epoch unix.
%
% See also FFT, IRF_SPECTROGRAM

ii = find(~isnan(data(:,1)));
if isempty(ii), error('time is NaN');
else ts = data(ii(1),1);
end
% Number of intervals must be computed from time
nint = 1;
ncomp = size(data,2) - 1; % first column is time
tcur = ts;


% Check if there is enough data
if( nint<1 )
  outF = [];
  outPxx = [];
  outspecrec = [];
  return
end

if nfft/2==fix(nfft/2), nf = nfft/2;
else nf = (nfft+1)/2;
end

specrec.f = sfreq*((1:nf) -1)'/nfft;
for jj=1:ncomp, specrec.p(jj) = {zeros(nint,nf)}; end
specrec.t = zeros(nint,1);

w = hanning(nfft);
wnorm = sum(w.^2)/nfft;	% normalization factor from windowing
nnorm = 2.0/nfft/sfreq/wnorm;

for jj=1:nint
  if(tseries)
    % Possibly FIXME (specrec.t to GenericTime once irf_spectrogram can handle TSeries)
    specrec.t(jj) = tcur.epochUnix + (nfft-1)/sfreq*.5;
    X = order_data(tlim(data, irf.tint(tcur, (nfft-1)/sfreq)), ...
      nfft, sfreq, tcur);
  else
    specrec.t(jj) = tcur + (nfft-1)/sfreq*.5;
    X = order_data(irf_tlim(data, tcur, tcur+(nfft-1)/sfreq), ...
      nfft, sfreq, tcur);
  end
  for comp=1:ncomp
    if( ~isempty(X) && any(~isnan(X(:,comp))) )
      ff = fft(detrend(X(:,comp)) .* w, nfft);
      pf = ff .*conj(ff) * nnorm;
      specrec.p{comp}(jj,:) = pf(1:nf);
    else
      specrec.p{comp}(jj,:) = NaN;
    end
  end
  % Next time interval (GenericTimeArray adds "double" as seconds)
  tcur = tcur + (1-overlap*.01)*(nfft-1)/sfreq;
end

if nargout==1, outspecrec = specrec;
elseif nargout==3,
  outspecrec = specrec.t;
  outPxx = specrec.p;
  outF = specrec.f;
elseif nargout==0,
  irf_spectrogram(specrec)
else
  error('irf_powerfft: unknown number of output parameters');
end

% Help function to clear datagaps
% We throw away intervals with less than 90% of data
function out = order_data(in, ndata, sfreq, ts)
  if(tseries)
    if isempty(in.data), out = []; return, end
    ncomp2 = size(in.data, 2);
    ind = round((in.time.epochUnix-ts.epochUnix)*sfreq+1);
    out = NaN(ndata, ncomp2);
    out(ind, :) = in.data;
  else
    if isempty(in), out = []; return, end
    ncomp2 = size(in,2) - 1; % Drop time column
    ind = round((in(:,1)-ts)*sfreq+1);
    out = NaN(ndata, ncomp2);
    out(ind, :) = in(:, 2:end); % Drop time column
  end
  indNaN = isnan(out);
  if( any(indNaN) )
    % if data has less than 90% return NaN, else replace with mean (excl. NaN).
    m = irf.nanmean(out, 1, 0.9);
    for col=1:ncomp2
      if(any(indNaN(:,col))), out(indNaN(:,col),col) = m(col); end;
    end
  end
end

end
