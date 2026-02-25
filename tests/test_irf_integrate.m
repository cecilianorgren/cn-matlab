% test_irf_integrate to see if it's necessary to correct for phase shift
%% Real data
cd /Users/Cecilia/Data/Cluster/20070417/;
load diE2;
diE2 = irf_abs(diE2);
t1 = irf_time('2007-04-17T15:32:50.666494131Z','utc>epoch');
t2 = irf_time('2007-04-17T15:32:50.802920627Z','utc>epoch');
tint = [t1 t2];
fmin = 40; fmax = 0;
acE = irf_filt(irf_tlim(diE2(:,[1 5]),tint),fmin,0);
tvec = acE(:,1);
% IRF integration
intE_a = irf_integrate(acE); 

% Matlab integration
dt = diff(acE(1:2,1));

intE_b = [acE(:,1) cumtrapz(acE(:,2))*dt];

% Plot and compare
h=irf_plot(1);
irf_plot(h,{[acE(:,1) acE(:,2)/max(acE(:,2))] [intE_a(:,1) intE_a(:,2)/max(intE_a(:,2))] [intE_b(:,1) intE_b(:,2)/max(intE_b(:,2))]},'comp')
irf_legend(h,{'E','irf\_integrate(E) (cumsum)','cumtrapz(E)'},[0.98, 0.95]);
irf_zoom(h,'x',tint)

%% Sinusoudal curve
nt = size(tvec,1);
dt = diff(acE(1:2,1));
step = 2;
phase = pi;
acE = [tvec(1:step:end) cos((1:step:nt)/pi-phase)'];
dacE = [tvec(1:step:end) sin((1:step:nt)/pi-phase)'];
% IRF integration
intE_a = irf_integrate(acE); 

% Matlab integration
intE_b = [tvec(1:step:end) cumtrapz(acE(:,2))*dt];

% Plot and compare
h=irf_plot(1);
irf_plot(h,{[acE(:,1) acE(:,2)/max(acE(:,2))],[dacE(:,1) dacE(:,2)/max(dacE(:,2))],[intE_a(:,1) intE_a(:,2)/max(intE_a(:,2))],[intE_b(:,1) intE_b(:,2)/max(intE_b(:,2))]},'comp')
irf_legend(h,{'E','int(E)','irf\_integrate(E) (cumsum)','cumtrapz(E)'},[0.98, 0.95]);
irf_zoom(h,'x',tint)