% While writing document
getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'a');
c_load('Atwo3')
c_load('wE3p32')
phase = c_phase(wE3p32(:,1),Atwo3);
tref = toepoch([2013 03 29 05 17 58.367]); % 10 min after st
h = irf_plot(phase);
irf_zoom(h,'x',tref + [-10 10]);
irf_pl_mark(tref + [-0.05 0.05])
c_pl_sc_orient(3,tref)

%%
getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'e')

%%
c_load('wE3p32')
c_load('wE3p34')
wE3p24 = irf_add(1,wE3p34,-1,wE3p32);
wE3p42 = irf_add(-1,wE3p34,1,wE3p32);
wE3p42plus = irf_add(1,wE3p34,1,wE3p32);

h = irf_plot({wE3p32,wE3p34,wE3p24,wE3p42,wE3p42plus});
ylabel(h(1),'wE3p32')
ylabel(h(2),'wE3p34')
ylabel(h(3),'wE3p34-wE3p32')
ylabel(h(4),'-wE3p34+wE3p32')
ylabel(h(5),'wE3p34+wE3p32')

%% See modulation frequencies
%specrec512 = irf_powerfft(wE3p32,2^9,c_efw_fsample(wE3p32));
%specrec1024 = irf_powerfft(wE3p32,2^10,c_efw_fsample(wE3p32));
%specrec2048 = irf_powerfft(wE3p32,2^11,c_efw_fsample(wE3p32));
%specrec4096 = irf_powerfft(wE3p32,2^12,c_efw_fsample(wE3p32));
specrec16384 = irf_powerfft(wE3p32,2^14,c_efw_fsample(wE3p32));
h = irf_spectrogram(specrec16384);
irf_zoom(h,'y',[0 2])
title('specrec16384')

