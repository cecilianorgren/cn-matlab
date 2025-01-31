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

mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

ic = 1; 
tint_load = irf.tint('2017-07-06T13:54:00.00Z/2017-07-06T13:55:00.000Z');
tint_load = irf.tint('2017-07-11T22:29:23.00Z/2017-07-11T22:30:23.000Z');
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint_load);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint_load);',ic);

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

tint_load = irf.tint('2017-07-11T22:29:23.00Z/2017-07-11T22:30:23.000Z');('2017-07-06T13:54:32.103127685Z/2017-07-06T13:54:32.315286621Z');
tint = irf.tint('2017-07-11T22:30:21.950Z/2017-07-11T22:30:21.980Z');
%tint = tint + [0.1 -0.1];
%tint = tint + [0.1 -0.1];
%%
tsEpar = gseE1par.tlim(tint);

y1 = tsEpar.data;
f_esw = @(t,tp2p) t.*exp(-(t./tp2p).^2);

t1 = tsEpar.time-tsEpar.time.start;
tcenter = median(t1);
dt = t1(2)-t1(1);

i1 = 1:numel(t1);

tp2p = dt*3; % sampling steps
t2 = -3*tp2p:dt:3*tp2p;
i2 = 1:numel(t2);

y2 = f_esw(t1-tcenter,tp2p);


z = xcorr(y1,y2);


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