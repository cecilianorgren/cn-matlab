S = @(t,f) cos(2*pi*f*t);

f = 2;
T = 1/f;
nT = 4;
nt = 2000;
t = linspace(0,nT*T,nt);

fs = 10;
ts1 = 0;
ts2 = nT*T;
ts = ts1:1/fs:ts2;
plot(t,S(t,f),ts,S(ts,f),'-')