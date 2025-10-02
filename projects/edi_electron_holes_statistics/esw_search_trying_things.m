%% Convolutions: 2 square functions

x = linspace(-1,1,100);
x1 = -1:0.02:1;
x2 = -0.1:0.02:0.1;

y1 = heaviside(x1-(0.89)).*(1-heaviside(x1-(0.99)));
y2 = heaviside(x2-(-0.09)).*(1-heaviside(x2-(0.07)));
z=conv(y1,y2);

if 0
hca = subplot(2,1,1);
plot(hca,1:numel(x),y1,1:numel(x),y2);
hca = subplot(2,1,2);
plot(hca,1:numel(z),z);
end

%

nRows = 4;
nCols = 1;
ip = 0;
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(iRow,iCol) = subplot(nRows,nCols,ip);
  end
end
isub = 1;

i1 = 1:numel(y1);
i2 = 1:numel(y2);
iz = 1:numel(z);

nz = numel(iz);
n1 = numel(i1);
n2 = numel(i2);

iref = 0;
for iref = 0:numel(z)-1
  %iref = 100;
  
  isub = 1;

  hca = h(isub); isub = isub + 1;
  plot(hca,i1,y1,i2,y2,'LineWidth',2); % ESW shape function
  
  hca = h(isub); isub = isub + 1;
  plot(hca,i1+n2,y1,'LineWidth',2); % data and ESW shape function

  hca = h(isub); isub = isub + 1;
  plot(hca,i2+iref,y2,'LineWidth',2); % data and ESW shape function
 
  
  hca = h(isub); isub = isub + 1;
  plot(hca,iz+n2/2,z,iref+1+n2/2,z(iref+1),'*','LineWidth',2);
  
  
  hlinks = linkprop(h(2:4),{'XLim'});
  h(2).XLim = [0 nz+n2];
  


  pause(0.0001)
end

%% Convolutions: Gaussian and square function

x = -10:.1:10;

y1 = exp(-x.^2);
y2 = heaviside(x-(-1)).*heaviside(x-(1));
y2 = heaviside(x-(-1)).*(1-heaviside(x-(1)));


z=conv(y1,y2);
hca = subplot(2,1,1);
plot(hca,x,y1,x,y2);
hca = subplot(2,1,2);
plot(hca,z);

%% Convolutions: 2 gaussian and gaussian

x = -10:.1:10;

y1 = exp(-(x+3).^2) + exp(-(x-3).^2);
y2 = exp(-x.^2);


z = conv(y1,y2);

hca = subplot(2,1,1);
plot(hca,x,y1,x,y2);
hca = subplot(2,1,2);
plot(hca,z);

%% Convolutions: Parallel electric field and dipolar function (deriv. of gaussian)

%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

ic = 1:4; 
tint_load = irf.tint('2017-07-06T13:54:00.00Z/2017-07-06T13:55:00.000Z');
tint_load = irf.tint('2017-07-11T22:29:23.00Z/2017-07-11T22:30:23.000Z');
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint_load);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint_load);',ic);

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

tint_load = irf.tint('2017-07-11T22:29:23.00Z/2017-07-11T22:30:23.000Z');('2017-07-06T13:54:32.103127685Z/2017-07-06T13:54:32.315286621Z');
tint = irf.tint('2017-07-11T22:30:21.950Z/2017-07-11T22:30:21.980Z');
%tint = tint + [0.1 -0.1];
tint = tint + 10*[-1 1];
%%
tsEpar = gseE1par.tlim(tint);

y1 = tsEpar.data;
f_esw = @(t,tp2p) t.*exp(-(t./tp2p).^2);

t1 = tsEpar.time-tsEpar.time.start;
tcenter = median(t1);
dt = t1(2)-t1(1);

i1 = 1:numel(t1);

tp2p = dt*3; % sampling steps
t2 = (-3*tp2p:dt:3*tp2p)';
i2 = 1:numel(t2);

y2 = f_esw(t2,tp2p);

nfitfun = numel(y2);


z = xcorr(y1,y2);
irem = [1:nfitfun-1 numel(c)-(nfitfun-1)+1:numel(c)];
z(irem) = [];
%%

[z_peaks,inds] = findpeaks(z,'MinPeakWidth',2);




c(irem) = [];
lags(irem) = [];
imax = lags(c==max(c));



tpeaks = inds - (numel(i2)+1)/2;
tsTesw = irf.ts_scalar(t1(tpeaks),ones(numel(tpeaks,1)));


i1 = 1:numel(y1);
i2 = 1:numel(y2);
iz = 1:numel(z);

nz = numel(iz);
n1 = numel(i1);
n2 = numel(i2);

nRows = 3;
nCols = 1;
ip = 0;
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(ip) = subplot(nRows,nCols,ip);
  end
end
isub = 1;

iref = 0;
for iref = 0:nz-1
  %iref = 100;
  
  isub = 1;

  hca = h(isub); isub = isub + 1;
  plotyy(hca,i1,y1,i2,y2); % ESW shape function
  
  %hca = h(isub); isub = isub + 1;
  %plot(hca,i1+n2,y1,'LineWidth',2); % data and ESW shape function

  %hca = h(isub); isub = isub + 1;
  %plot(hca,i2+iref,y2,'LineWidth',2); % data and ESW shape function
  
  hca = h(isub); isub = isub + 1;
  plotyy(hca,i1+n1,y1,i1+iref,y2); % data and ESW shape function


  hca = h(isub); isub = isub + 1;
  plot(hca,iz+n1,z,iref+1+n2/2,z(iref+1),'*','LineWidth',2);
  
  
  hlinks = linkprop(h(2:3),{'XLim'});
  h(2).XLim = [0 nz+n2];
  


  pause(0.0001)
end

%%
for iref = 0:numel(z)-1
  %iref = 100;
  
  isub = 1;

  hca = h(isub); isub = isub + 1;
  plot(hca,i2,y2); % ESW shape function
  
  hca = h(isub); isub = isub + 1;
  plotyy(hca,i2+iref,y2,i1+(numel(i2)-1),y1); % data and ESW shape function
  
  hca = h(isub); isub = isub + 1;
  plot(hca,numel(i2)+(1:numel(z)),z,iref+numel(i2),z(iref+1),'*');
  
  
  hlinks = linkprop(h(2:3),{'XLim'});
  h(2).XLim = [0 numel(z)+2*numel(y2)];
  


  pause(0.01)
end

%%



x = 0:1000;
x1 = -10:10;
y1 = x1.*exp(-((x1)/5).^2);
x2 = -20:20;
y2 = (x2-4).*exp(-((x2-4)/5).^2);

i1 = 1:numel(x1);
i2 = 1:numel(x2);

[c,lags] = xcorr(y2,y1);
imax = lags(c==max(c));

%plot(x1+numel(x1)+1-imax,y1,x2,y2,lags,c12)
plot(i2,y2,i1+imax,y1,imax+(numel(i1)+1)/2,0,'sk')

%%
x = 0:1000;
x1 = -10:10;
y1 = x1.*exp(-((x1)/5).^2);
x2 = -20:20;
y2 = (x2-4).*exp(-((x2-4)/5).^2);

nfitfun = numel(x1);

i1 = 1:numel(x1);
i2 = 1:numel(x2);

[c,lags] = xcorr(y2,y1);
irem = [1:nfitfun-1 numel(c)-(nfitfun-1)+1:numel(c)];
c(irem) = [];
lags(irem) = [];
imax = lags(c==max(c));

%plot(x1+numel(x1)+1-imax,y1,x2,y2,lags,c12)
plot(i2,y2,i1+imax,y1,imax,0,'o')





%% Find peaks in int(Epar), then make fit based on peak locations
tsEpar = gseE1par.tlim(tint(1)+[0  10]);
timeline = tsEpar.time;

tsEpar_filt = tsEpar.filt(10,0,[],5);

y1 = tsEpar.data;
f_esw = @(t,tp2p) t.*exp(-(t./tp2p).^2);


dt = timeline(2) - timeline(1);
t_min_esw = 0.001;
i_min_esw = t_min_esw/dt;
E_thres = 2; % mV/m

intE = irf_integrate(tsEpar);
intEfilt =  intE.filt(10,0,[],5);
negintEfilt = 1*intEfilt; % do both on negative and positive peaks... but not on abs(), becuse then it just becomes really messy
%[peaks,inds] = findpeaks(negintEfilt.data,'MaxPeakWidth',0.5*i_min_esw);
MinPeakWidth = 0.5*i_min_esw;
[peaks,inds, width, prom] = findpeaks(negintEfilt.data,'MinPeakWidth',MinPeakWidth,'MinPeakDistance',MinPeakWidth*4);
minPeakProminence = prctile(prom,75);
[peaks, inds, width, prom] = findpeaks(negintEfilt.data,'MinPeakWidth',MinPeakWidth,'MinPeakDistance',MinPeakWidth*4,'MinPeakProminence',minPeakProminence);

all_esw = struct([]);

% Sort out intervals with E_max < E_thres
iskip = [];
for ipeak = 1:numel(inds)
  if mod(ipeak,100) == 0 ; disp(sprintf('%g/%g',ipeak,numel(inds))); end
  all_esw(ipeak).peak.time = tsEpar_filt.time(inds(ipeak));
  all_esw(ipeak).peak.width = width(ipeak);
  all_esw(ipeak).peak.prom = prom(ipeak);
  Etint = all_esw(ipeak).peak.time + 0.5*all_esw(ipeak).peak.width*dt*[-1 1];
  Etint;
  Etmp = tsEpar.tlim(Etint).data;
  if max(abs(Etmp)) < E_thres
    %all_esw(ipeak) = [];
    iskip(end+1) = ipeak;
  end
end
all_esw(iskip) = [];
nESW = numel(all_esw);
%%

% t_width = t_min_esw*5;
skippeak = [];
for ipeak = 1:nESW%numel(inds)
  if mod(ipeak,100) == 0 ; disp(sprintf('%g/%g',ipeak,nESW)); end
  max_t_width = dt*width(ipeak)*3;
  ttmp = timeline(inds(ipeak));
  Etmp = tsEpar.tlim(ttmp + max_t_width*0.5*[-1 1]);
  [fitparam,tsESW,exitflag] = fit_esw(Etmp);
  if exitflag == 0
    skippeak(end+1) = ipeak;
  end
  all_esw(ipeak).time = Etmp.time.start + fitparam(3);
  all_esw(ipeak).tp2p = fitparam(2);
  all_esw(ipeak).phi0 = fitparam(3);
  all_esw(ipeak).offset = fitparam(4);
  all_esw(ipeak).tsESW = tsESW;
  all_esw(ipeak).peak.time = tsEpar.time(inds(ipeak));
  all_esw(ipeak).peak.width = width(ipeak);
  all_esw(ipeak).peak.prom = prom(ipeak);
  %pause(1)
  

end
all_esw(skippeak) = [];
%
tsESW_findpeaks = irf.ts_scalar(tsEpar.time(inds),peaks);
%tsESW_fit = irf.ts_scalar([all_esw.time],peaks);

%%
h = irf_plot(4);

hca = irf_panel('Epar');
irf_plot(hca,{tsEpar})
hca.YLabel.String = 'E_{||} (mV/m)';

hold(hca,'on')
for iESW = 1:numel(all_esw)
  irf_plot(hca,all_esw(iESW).tsESW,'r')
end
hold(hca,'off')

hca = irf_panel('int Epar');

irf_plot(hca,{intE})
hca.YLabel.String = '\int E dt (s*mV/m)';

hca = irf_panel('int Epar filt');
irf_plot(hca,{intEfilt})
%hold(hca,'on')
%irf_plot(hca,{tsESW},'*')
%hold(hca,'off')
hca.YLabel.String = '\int E_{filt}dt (s*mV/m)';

hca = irf_panel('neg int Epar filt');
irf_plot(hca,{negintEfilt})
hold(hca,'on')
irf_plot(hca,{tsESW_findpeaks},'*')
irf_plot(hca,{tsESW_fit},'s')
hold(hca,'off')
hca.YLabel.String = '-\int E_{filt}dt (s*mV/m)';



%% Find peaks in int(Epar), then make fit based on peak locations, 2
tsEpar = gseE1par.tlim(tint(1)+[10  30]);
timeline = tsEpar.time;

tsEpar_filt = tsEpar.filt(10,0,[],5);

y1 = tsEpar.data;
f_esw = @(t,tp2p) t.*exp(-(t./tp2p).^2);


dt = timeline(2) - timeline(1);
t_min_esw = 0.001;
i_min_esw = t_min_esw/dt;
E_thres = 2; % mV/m

intE = irf_integrate(tsEpar);
intEfilt =  intE.filt(10,0,[],5);
negintEfilt = 1*intEfilt; % do both on negative and positive peaks... but not on abs(), becuse then it just becomes really messy
%[peaks,inds] = findpeaks(negintEfilt.data,'MaxPeakWidth',0.5*i_min_esw);
MinPeakWidth = 0.5*i_min_esw;
[peaks,inds, width, prom] = findpeaks(negintEfilt.data,'MinPeakWidth',MinPeakWidth,'MinPeakDistance',MinPeakWidth*4);
minPeakProminence = prctile(prom,75);
[peaks, inds, width, prom] = findpeaks(negintEfilt.data,'MinPeakWidth',MinPeakWidth,'MinPeakDistance',MinPeakWidth*4,'MinPeakProminence',minPeakProminence);

all_esw = struct([]);

% Sort out intervals with E_max < E_thres
iskip = [];
for ipeak = 1:numel(inds)
  if mod(ipeak,100) == 0 ; disp(sprintf('%g/%g',ipeak,numel(inds))); end
  all_esw(ipeak).peak.time = tsEpar_filt.time(inds(ipeak));
  all_esw(ipeak).peak.width = width(ipeak);
  all_esw(ipeak).peak.prom = prom(ipeak);
  Etint = all_esw(ipeak).peak.time + 0.5*all_esw(ipeak).peak.width*dt*[-1 1];
  Etint;
  Etmp = tsEpar.tlim(Etint).data;
  if max(abs(Etmp)) < E_thres
    %all_esw(ipeak) = [];
    iskip(end+1) = ipeak;
  end
end
all_esw(iskip) = [];
nESW = numel(all_esw);
%%
% t_width = t_min_esw*5;
for ipeak = 1:nESW%numel(inds)
  if mod(ipeak,100) == 0 ; disp(sprintf('%g/%g',ipeak,nESW)); end
  max_t_width = dt*width(ipeak)*3;
  ttmp = timeline(inds(ipeak));
  Etmp = tsEpar.tlim(ttmp + max_t_width*0.5*[-1 1]);
  [fitparam,tsESW] = fit_esw(Etmp);
  all_esw(ipeak).time = Etmp.time.start + fitparam(3);
  all_esw(ipeak).tp2p = fitparam(2);
  all_esw(ipeak).phi0 = fitparam(3);
  all_esw(ipeak).offset = fitparam(4);
  all_esw(ipeak).tsESW = tsESW;
  all_esw(ipeak).peak.time = tsEpar.time(inds(ipeak));
  all_esw(ipeak).peak.width = width(ipeak);
  all_esw(ipeak).peak.prom = prom(ipeak);
  %pause(1)
  

end
%
tsESW_findpeaks = irf.ts_scalar(tsEpar.time(inds),peaks);
%tsESW_fit = irf.ts_scalar([all_esw.time],peaks);

