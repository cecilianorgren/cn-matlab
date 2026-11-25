% powerfft and TSpectrs

fftE = irf_powerfft(dslE3,256,8192);

newTime = EpochTT(fftE.t);

Ts = TSpectra(newTime,fftE.p,fftE.f,dslE3.representation);