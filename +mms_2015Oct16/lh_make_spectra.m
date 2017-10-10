tint_fft = irf.tint('2015-10-16T10:33:24.00Z',2); % magnetosphere-magnetosheath-magnetosphere
ic = 3;

c_eval('locE? = dslE?brst.tlim(tint_fft).dot(direction);',ic)
c_eval('fs = 1/(locE?.time(2)-locE?.time(1))',ic);
c_eval('pfftE? = irf_powerfft(locE?,locE?.length,fs,0.5);',ic)

c_eval('locB? = dmpaB?scm.tlim(tint_fft).dot(direction);',ic)
c_eval('fs = 1/(locB?.time(2)-locB?.time(1))',ic);
c_eval('pfftB? = irf_powerfft(locB?,locB?.length,fs,0.5);',ic)

fn = 100; % number of f intervals for smoothing
cn_smooth_fft(plotEfft,fn)

phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
v_factor = 0.5;
phiEscale = velocity*v_factor*1e3;
di_loc = 1;
    
c_eval('pfftE = pfftE?',ic)
c_eval('pfftB = pfftB?',ic)
hca = subplot(1,1,1);
lines = loglog(hca,pfftB.f*2*pi/velocity*di_loc,pfftB.p{1}*1e-9*phiBscale^2,...
                   pfftE.f*2*pi/velocity*di_loc,pfftE.p{1}./pfftE.f'*phiEscale^2);