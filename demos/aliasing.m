S = @(t,f,thsift) cos(2*pi*f*(t-tshift));

% real signal
f = 2;
T = 1/f;
nT = 10;
nt = 2000;
t = linspace(0,nT*T,nt);

fs = 1.2*f; % sampling frequency needs to be atleast 2*f
tshift = 0.2/f;
ts1 = 0;%0.1/f; % add a small time shift so that we do not shift the
ts2 = nT*T;
ts = ts1:1/fs:ts2;

plot(t,S(t,f,tshift),'-',ts,S(ts,f,tshift),'*-','linewidth',1.5)
legend({sprintf('real signal, f = %g',f),sprintf('sampled signal, fs = %g = %gf',fs,fs/f)},'location','northoutside')
title('Misinterpreting the signal if the sampling frequency is too low')