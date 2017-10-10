%% fft
fs=512;
t=1/fs:1/fs:1;
ft=linspace(1/fs,fs/8,fs);
f0=150;
B=1; % nT
theta1=2*pi*ft;
theta2=2*pi*f0;
ut=sin(theta1.*t)+sin(theta2*t);
subplot(3,1,1)
plot(t,ut)
%%
uw=fft(ut);
subplot(3,1,2)
plot(1:fs/2,uw(1:fs/2))
%%
uw2=fft(ut(end/2:end/2+63))
subplot(3,1,3)
plot(1:64,uw2)
%%
hamming_w=@(j,L)(0.54-0.46*cos(2*pi*j./(L-1)));
hann_w=@(j,L)(0.5-0.5*cos(2*pi*j'./(L-1)));
gaussian=@(j,L,K)(exp(-(j-0.5*(L-1)).^2./(2*K.^2)));
j=0:0.01:2*pi;
%%
t=-5:0.01:5;
f0=2*pi/10;
morl=@(t,f0)(exp(-(t).^2/2).*exp(-1i*f0*(t)));

