function out = cn_smooth_fft(specrec,nf)
% Make spectra smooth by averaging.
% First run irf_powerfft.
%nf = 100; % number of f intervals to be divided logarithmically
out = specrec;
f = specrec.f;
p = specrec.p{1};
f1 = f(1); f2 = f(end);
logf = log10(f);

fEdges = linspace(logf(2),logf(end),nf+1);
fCenter = fEdges(1:end-1)+0.5*(fEdges(2)-fEdges(1));

[N,binIdx] = histc(logf, fEdges);
meanp = nan(nf,1);
for ii = 1:nf
    pp = mean(p(find(binIdx==ii)));
    if isempty(pp)
        pp = NaN;
    end
    meanp(ii) = pp;        
end    
out.f = 10.^fCenter; 
out.p = {meanp};
