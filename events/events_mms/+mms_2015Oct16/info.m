%% The whole time interval
tint = irf.tint('2015-10-16T05:30:00.00Z/2015-10-16T16:00:00.00Z');

%% Cold ions
tint = irf.tint('2015-10-16T14:50:00.00Z/2015-10-16T14:57:00.00Z');
 
%% Cold ions subinterval with very sharp crossing
tint = irf.tint('2015-10-16T14:50:00.00Z/2015-10-16T14:57:00.00Z');

%% Candidate diffusion region
tint = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z');
tint = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z');
tDRCenter = irf.tint('2015-10-16T10:33:28.00Z',0.00001); tDRCenter = tDRCenter(1);

%% Candidate diffusion region, entire out-in-out crossing
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');

%% Candidate diffusion region, entire out-in-out-in(?) crossing
tint = irf.tint('2015-10-16T10:33:00.00Z/2015-10-16T10:36:00.00Z');
%% Roy Torbert
tint = irf.tint('2015-10-16T13:00:00',8*60);

%% Roy Torbert: really deep spacecraft potential
tint = irf.tint('2015-10-16T11:25:00',60*7);


irf_plot({ne1,tepar1,teper1})

B = 20, n = 15, Te = 100, Ti = 300;

Le = 1.37
Li = 59
rhoe = 3
rhop = 250