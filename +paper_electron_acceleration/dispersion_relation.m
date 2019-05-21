ic = 1;

%% Set tint
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
tint_fig = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:54:45.00Z');
tint_lobe = irf.tint('2017-07-06T00:54:07.00Z/2017-07-06T00:54:08.00Z');
tint_sheet = irf.tint('2017-07-06T00:54:18.70Z/2017-07-06T00:54:19.50Z');
tint_sep = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:16.50Z');
tint_phi = irf.tint('2017-07-06T00:54:14.00Z/2017-07-06T00:54:15.50Z');
tint_fred = tint_phi+[-1 1];
tint_disprel_slow = irf.tint('2017-07-06T00:54:14.20Z/2017-07-06T00:54:14.40Z');
tint_disprel_fast = irf.tint('2017-07-06T00:54:15.20Z/2017-07-06T00:54:15.40Z');
times_fred = irf_time('2017-07-06T00:54:14.20Z','utc>epochtt') + [0 0.1 0.2 0.3];
  
%% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
localuser = datastore('local','user');
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS/']);
db_info = datastore('mms_db');   
pathLocalUser = ['/Users/' localuser '/'];

%% Load data
units = irf_units;
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Remove all one-count "noise"
%c_eval('iPDist?_no1c.data(iPDist?_no1c.data<iPDistErr?.data*1.1) = 0;',ic)
c_eval('iPDist?_no1c = iPDist?;',ic)
c_eval('iPDist?_no1c.data(iPDist?_no1c.data<iPDistErr?.data*1.1) = 0;',ic)
c_eval('ePDist?_no1c = ePDist?;',ic)
c_eval('ePDist?_no1c.data(ePDist?_no1c.data<ePDistErr?.data*1.1) = 0;',ic)

%% Choose spacecraft
ic = 1;
c_eval('eDist = ePDist?.tlim(tint_fred);',ic)
c_eval('iDist = iPDist?.tlim(tint);',ic)

c_eval('eDist_no1c = ePDist?_no1c.tlim(tint_fred);',ic)
c_eval('iDist_no1c = iPDist?_no1c.tlim(tint);',ic)


% Smooth electron distribution, sliding average
eDist_sm3 = eDist;
ed1 = eDist.resample(eDist.time(1:3:end));
ed2 = eDist.resample(eDist.time(2:3:end));
ed3 = eDist.resample(eDist.time(3:3:end));
eDist_sm3.data(1:3:end,:) = ed1.data(:,:);
eDist_sm3.data(2:3:end,:) = ed2.data(:,:);
eDist_sm3.data(3:3:end,:) = ed3.data(:,:);

eDist_no1c_sm3 = eDist_no1c;
ed1 = eDist.resample(eDist_no1c.time(1:3:end));
ed2 = eDist.resample(eDist_no1c.time(2:3:end));
ed3 = eDist.resample(eDist_no1c.time(3:3:end));
eDist_no1c_sm3.data(1:3:end,:) = ed1.data(:,:);
eDist_no1c_sm3.data(2:3:end,:) = ed2.data(:,:);
eDist_no1c_sm3.data(3:3:end,:) = ed3.data(:,:);

%% Make reduced distribution, remove background
ZeroNaN = 0;
tic; [eDist1_bgremoved, eDist1_bg, ephoto_scale] = mms.remove_edist_background(ePDist1, 'tint', tint_fred, 'ZeroNaN', ZeroNaN); toc
eDist = eDist1_bgremoved;  
eDist_bgremoved = eDist1_bgremoved;  
eDist_original = ePDist1;

%% Make reduced distribution, reduce distributon
scpot = scPot1.resample(eDist);
lowerelim = scpot*1.3;
par_dir = dmpaB1.resample(eDist).norm;
% energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
% vg = [-energies(end:-1:1) 0 energies];
% vg(abs(vg)>80000) = [];
vg_e = -50000:1000:50000;
vg_i = -2500:20:2500;
tic; ef1D = eDist.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scpot); toc % reduced distribution along B
tic; ef1D_no1c = eDist_no1c.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scpot); toc % reduced distribution along B
tic; ef1D_no1c_sm3 = eDist_no1c_sm3.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scpot); toc % reduced distribution along B
tic; ef1D_sm3 = eDist_sm3.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scpot); toc % reduced distribution along B
tic; if1D_long = iPDist1.reduce('1D',dmpaB1.resample(iPDist1).norm,'lowerelim',lowerelim,'vg',vg_i); toc % reduced distribution along B
tic; if1D_no1c_long = iPDist1_no1c.reduce('1D',dmpaB1.resample(iPDist1).norm,'lowerelim',lowerelim,'vg',vg_i); toc % reduced distribution along B
tic; if1D = iDist.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_i); toc % reduced distribution along B
tic; if1D_no1c = iDist_no1c.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_i); toc % reduced distribution along B


if1D_no1c_long_nan = if1D_no1c_long;
if1D_no1c_long_nan.data(if1D_no1c_long_nan.data<1e-4) = NaN; 

tic; ef1D_no1c_long = ePDist1_no1c.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scPot1.resample(ePDist1_no1c)); toc % reduced distribution along B
ef1D_no1c_long_nan = ef1D_no1c_long;
ef1D_no1c_long_nan.data(ef1D_no1c_long_nan.data<1e-6) = NaN; 

tic; ef1D_long = ePDist1.reduce('1D',par_dir,'lowerelim',lowerelim,'vg',vg_e,'scpot',scPot1.resample(ePDist1_no1c)); toc % reduced distribution along B
ef1D_long_nan = ef1D_long;
ef1D_long_nan.data(ef1D_long_nan.data<1e-6) = NaN; 

% smooth
eDist_no1c_sm3 = eDist_no1c;
ed1 = eDist.resample(eDist_no1c.time(1:3:end));
ed2 = eDist.resample(eDist_no1c.time(2:3:end));
ed3 = eDist.resample(eDist_no1c.time(3:3:end));
eDist_no1c_sm3.data(1:3:end,:) = ed1.data(:,:);
eDist_no1c_sm3.data(2:3:end,:) = ed2.data(:,:);
eDist_no1c_sm3.data(3:3:end,:) = ed3.data(:,:);

ef1D_no1c_long_sm3 = ef1D_no1c_long.smooth(3);
ef1D_no1c_long_sm5 = ef1D_no1c_long.smooth(5);

ef1D_no1c_long_sm3_nan = ef1D_no1c_long_sm3;
ef1D_no1c_long_sm3_nan.data(ef1D_no1c_long_sm3.data<1e-6) = NaN; 

ef1D_no1c_long_sm5_nan = ef1D_no1c_long_sm5;
ef1D_no1c_long_sm5_nan.data(ef1D_no1c_long_sm5.data<1e-6) = NaN; 

%% Simple plot, compare long ions to electrons, for dispersion analysis purpose
h = irf_plot(6);
isub = 1;

hca = irf_panel('iPDist1');
irf_spectrogram(hca,iPDist1.deflux.omni.specrec)
hca.YScale = 'log';

hca = irf_panel('if1D_long');
irf_spectrogram(hca,if1D_long.specrec)


hca = irf_panel('iPDist1_no1c');
irf_spectrogram(hca,iPDist1_no1c.deflux.omni.specrec)
hca.YScale = 'log';

hca = irf_panel('if1D_no1c_long');
irf_spectrogram(hca,if1D_no1c_long.specrec)

%hca = irf_panel('if1D_no1c_long_nan');
%irf_spectrogram(hca,if1D_no1c_long_nan.specrec)

hca = irf_panel('ePDist1');
irf_spectrogram(hca,ePDist1.deflux.omni.specrec)
hca.YScale = 'log';

hca = irf_panel('ef1D_long');
irf_spectrogram(hca,ef1D_long.specrec)
%hca.CLim = [-6 -3];
hlink = linkprop(h,{'XLim'});

%% Simple plot to compare, long
h = irf_plot(8);
isub = 1;

hca = irf_panel('ePDist1');
irf_spectrogram(hca,ePDist1.omni.specrec)
%hca.CLim = [-6 -3];
hca.YScale = 'log';

hca = irf_panel('ePDist1.dpflux');
irf_spectrogram(hca,ePDist1.dpflux.omni.specrec)
%hca.CLim = [-6 -3];
hca.YScale = 'log';

hca = irf_panel('ef1D');
irf_spectrogram(hca,ef1D_long.specrec)
hca.CLim = [-6 -3];

hca = irf_panel('ef1D_nan');
irf_spectrogram(hca,ef1D_long_nan.specrec)
hca.CLim = [-6 -3];

hca = irf_panel('ePDist1_no1c');
irf_spectrogram(hca,ePDist1_no1c.omni.specrec)
%hca.CLim = [-6 -3];
hca.YScale = 'log';

hca = irf_panel('ePDist1_no1c.dpflux');
irf_spectrogram(hca,ePDist1_no1c.dpflux.omni.specrec)
%hca.CLim = [-6 -3];
hca.YScale = 'log';

hca = irf_panel('ef1D_no1c');
irf_spectrogram(hca,ef1D_no1c_long.specrec)
hca.CLim = [-6 -3];

hca = irf_panel('ef1D_no1c_nan');
irf_spectrogram(hca,ef1D_no1c_long_nan.specrec)
hca.CLim = [-6 -3];

if 0
  hca = irf_panel('ef1D_no1c_sm3');
  irf_spectrogram(hca,ef1D_no1c_long_sm3.specrec)
  hca.CLim = [-6 -3];
end
if 0
  hca = irf_panel('ef1D_no1c_sm3_nan');
  irf_spectrogram(hca,ef1D_no1c_long_sm3_nan.specrec)
  hca.CLim = [-6 -3];
end
if 0
  hca = irf_panel('ef1D_no1c_sm5');
  irf_spectrogram(hca,ef1D_no1c_long_sm5.specrec)
  hca.CLim = [-6 -3];
end
if 0
  hca = irf_panel('ef1D_no1c_sm5_nan');
  irf_spectrogram(hca,ef1D_no1c_long_sm5_nan.specrec)
  hca.CLim = [-6 -3];
end
%hlink = linkprop(h,{'XLim','YLim','CLim'});
hlink = linkprop(h,{'XLim'});

%% Simple plot to compare, short
h = irf_plot(6);
isub = 1;

hca = irf_panel('if1D');
irf_spectrogram(hca,if1D.specrec)

hca = irf_panel('if1D_no1c');
irf_spectrogram(hca,if1D_no1c.specrec)

hca = irf_panel('ef1D');
irf_spectrogram(hca,ef1D.specrec)

hca = irf_panel('ef1D_no1c');
irf_spectrogram(hca,ef1D_no1c.specrec)

hca = irf_panel('ef1D_sm3');
irf_spectrogram(hca,ef1D_sm3.specrec)

hca = irf_panel('ef1D_no1c_sm3');
irf_spectrogram(hca,ef1D_no1c_sm3.specrec)

%% Plot, electron distribution and fit

% fit
units = irf_units;% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps0 = 8.854e-12;
kB = 1.381e-23;
c = 299792458;

% 1D reduced Maxwellian distributions
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2).*exp(-(v-vd).^2/vt.^2);

% Measured
shift_slow = 0;
shift_fast = 0.05 + 0.00*[1 -1] ;
ffe_slow = ef1D.tlim(tint_disprel_slow + shift_slow).data;
ffe_fast = ef1D.tlim(tint_disprel_fast + shift_fast).data;
vve = ef1D.tlim(tint_disprel_slow).depend{1}(1,:)*1e-3;

h = setup_subplots(3,1);
isub = 1;

if 1 % slow
  hca = h(isub); isub = isub + 1;
  plot(hca,vve,mean(ffe_slow,1),'-')
  hca.XTick = -40:20:40;
  hca.XGrid = 'on';
  %hca.XLim = 
  % fit
  % Plasma properties for the different species
  [n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('slow_e');
  nsp = numel(n); % number of species
  ntot = sum(n); wptot = sqrt(ntot.*q(1).^2./(m(1)*eps0)); % Hz
  hold(hca,'on')
  plot(hca,vve,f(vve*1e6,n(1),vt(1),vd(1))+f(vve*1e6,n(2),vt(2),vd(2)))  
  for ipop = 1:2
    plot(hca,vve,f(vve*1e6,n(ipop),vt(ipop),vd(ipop)))
  end
  hold(hca,'off')
  hca.XTick = -40:20:40;
  hca.XLabel.String = [];
  hca.XGrid = 'on';
end

if 1 % fast
  hca = h(isub); isub = isub + 1;
  plot(hca,vve,mean(ffe_fast,1),'-')
  hca.XTick = -40:20:40;
  hca.XGrid = 'on';
  %hca.XLim = 
  % fit
  % Plasma properties for the different species
  [n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('fast_e');
  nsp = numel(n); % number of species
  ntot = sum(n); wptot = sqrt(ntot.*q(1).^2./(m(1)*eps0)); % Hz
  hold(hca,'on')
  plot(hca,vve,f(vve*1e6,n(1),vt(1),vd(1))+f(vve*1e6,n(2),vt(2),vd(2)))  
  for ipop = 1:2
    plot(hca,vve,f(vve*1e6,n(ipop),vt(ipop),vd(ipop)))
  end
  hold(hca,'off')
  hca.XTick = -40:20:40;
  hca.XLabel.String = [];
  hca.XGrid = 'on';
end

%% Previous
if 0 % timseries
  hca = h(isub); isub = isub + 1;
  tsf = ef1D;
  tsf.data(tsf.data<1e-6) = NaN;
  irf_spectrogram(hca,tsf.specrec)
  hca.YLim = 55*1e3*[-1 1];
  irf_pl_mark(hca,tint_disprel_slow+shift_slow,'r')
  irf_pl_mark(hca,tint_disprel_fast+shift_fast,'g')
  irf_timeaxis(hca)
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,vve,mean(ffe_slow,1),vve,mean(ffe_fast,1))
  hca.XTick = -40:20:40;
  hca.XGrid = 'on';
  %hca.XLim = 
end

hca = h(isub); isub = isub + 1;
if 1 % fit
  plot(hca,vve,f(vve*1e6,n(1),vt(1),vd(1))+f(vve*1e6,n(2),vt(2),vd(2)))
  hold(hca,'on')
  for ipop = 1:2
    plot(hca,vve,f(vve*1e6,n(ipop),vt(ipop),vd(ipop)))
  end
  hold(hca,'off')
  hca.XTick = -40:20:40;
  hca.XLabel.String = [];
  hca.XGrid = 'on';
end

if 0
  hca = h(isub); isub = isub + 1;
  ffi_1 = if1D.tlim(tint_disprel_slow+0.1*[-1 1]).data;
  ffi_2 = if1D.tlim(tint_disprel_slow+0.2*[-1 1]).data;
  vvi = if1D.tlim(tint_disprel_slow).depend{1}(1,:);
  plot(hca,vvi,mean(ffi_1,1),vvi,mean(ffi_2,1))
end

%hca = h(isub); isub = isub + 1;
%plot(hca,vve,mean(ffe_slow,1))

if 0
  hca = h(isub); isub = isub + 1;
  ffi_1 = if1D.tlim(tint_disprel_fast+0.1*[-1 1]).data;
  ffi_2 = if1D.tlim(tint_disprel_fast+0.2*[-1 1]).data;
  vvi = if1D.tlim(tint_disprel_fast).depend{1}(1,:);
  plot(hca,vvi,mean(ffi_1,1),vvi,mean(ffi_2,1))
end