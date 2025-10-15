%% Local coordinate system
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn_edr = [L;M;N];

L_vi = -[-0.8906    0.4548    0.0045];
M_vi = [ 0.4539    0.8893   -0.0559];
N_vi = -[-0.0294   -0.0477   -0.9984];
lmn_vi = [L_vi; M_vi; N_vi];


L_gse = [1 0 0];
M_gse = [0 1 0];
N_gse = [0 0 1];
lmn_gse = [L_gse; M_gse; N_gse];

L_gse = [1 0 -.2]; L_gse = L_gse/norm(L_gse);
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L_gse)); M_gse = M_gse/norm(M_gse);
%N_gse = [0 0 1];
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

lmn = lmn_vi;
lmn = lmn_gse;
%lmn = lmn_edr;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);


c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)

%% Prepare PDist iwht movmean
nMovMean = 7;
c_eval('pdist_cleaned = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)

%% Make PDist objects
% Test if the 360 deg spectrogram works
PD_orig = pdist_cleaned;
PD_shift_eye = PD_orig.shift_many([-170 0 0],200,[1 0 0; 0 1 0; 0 0 1],ic);
PD_shift_lmn = PD_orig.shift_many([-170 0 0],200,lmn,ic);
PD_shift_nlm = PD_orig.shift_many([-170 0 0],200,circshift(lmn,1,1),ic);
PD_shift_mnl = PD_orig.shift_many([-170 0 0],200,circshift(lmn,2,1),ic);
PD_shift_mnl_large = PD_orig.shift_many([-2000 0 0],200,circshift(lmn,2,1),ic);
PD_noshift_mnl = PD_orig.shift_many([0 0 0],200,circshift(lmn,2,1),ic);

dir_pitch = irf.ts_vec_xyz(iPDist3.time,repmat(N,iPDist3.length,1));
iPitch_orig  = PD_orig.pitchangles(dir_pitch,12);
iPitch_shift_eye = PD_shift_eye.pitchangles(dir_pitch,12);
iPitch_shift_lmn = PD_shift_lmn.pitchangles(dir_pitch,12);
iPitch_shift_nlm = PD_shift_nlm.pitchangles(dir_pitch,12);
iPitch_shift_mnl = PD_shift_mnl.pitchangles(dir_pitch,12);
iPitch_shift_mnl_large = PD_shift_mnl_large.pitchangles(dir_pitch,12);
iPitch_noshift_mnl = PD_noshift_mnl.pitchangles(dir_pitch,12);

iAzim_orig  = PD_orig.azimuthangle;
iAzim_shift_eye = PD_shift_eye.azimuthangle;
iAzim_shift_lmn = PD_shift_lmn.azimuthangle;
iAzim_shift_nlm = PD_shift_nlm.azimuthangle;
iAzim_shift_mnl = PD_shift_mnl.azimuthangle;
iAzim_shift_mnl_large = PD_shift_mnl_large.azimuthangle;
iAzim_noshift_mnl = PD_noshift_mnl.azimuthangle;

%% Plot
elim = [1500 Inf];

h = irf_plot(8);


if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1
  hca = irf_panel('azim original');
  specrec = iAzim_orig.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted eye');
  specrec = iAzim_shift_eye.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted lmn');
  specrec = iAzim_shift_lmn.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted nlm');
  specrec = iAzim_shift_nlm.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted mnl');
  specrec = iAzim_shift_mnl.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted mnl large');
  specrec = iAzim_shift_mnl_large.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim not shifted mnl');
  specrec = iAzim_noshift_mnl.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end

colormap(pic_colors('candy4'))
irf_plot_axis_align(h)

%% Macroparticles
units = irf_units;
PD = pdist_cleaned.elim([500 Inf]);
%PD = ;

nt = PD.length;
Vsc = scPot3.resample(PD);
nMP = 5000;
MP = PD.macroparticles('ntot',nMP,'scpot',Vsc);
vL_shift = -170*1;
%vL_shift = -999;
E_edges = [PD.ancillary.energy(1,1) - PD.ancillary.delta_energy_minus(1,1), PD.ancillary.energy(1,:) + PD.ancillary.delta_energy_plus(1,:)];
nEnergy = numel(E_edges)-1;
azimuth_edges = -180:5:180;
nAzimuth = numel(azimuth_edges)-1;

dn_tot_all_ = zeros(nt,nAzimuth);
dn_tot_all = zeros(nt,nEnergy,nAzimuth);
df_tot_all = zeros(nt,nEnergy,nAzimuth);
dv_tot_all = zeros(nt,nEnergy,nAzimuth);
for it = 1:nt
  elow = tsElow.resample(PD(it).time).data;
  dv = MP(it).dv;
  df = MP(it).df;
  dn = df.*dv;
  vx = MP(it).vx;
  vy = MP(it).vy;
  vz = MP(it).vz;
  vL = [vx vy vz]*L' - vL_shift;
  vM = [vx vy vz]*M';
  vN = [vx vy vz]*N';
  theta_NL = atan2d(vN,vL);
  v2 = (vL.^2 + vM.^2 + vN.^2)*1e6; % (m/s)^2
  vLN = sqrt(vL.^2 + vN.^2); % (km/s)^2
  E = units.mp*v2/2/units.eV; % eV
  E(E<elow) = 0;
  %dn(vM<0) = 0;
  %dn(vLN<1500) = 0;  
  
  [dn_tot_ edges_ mid_ loc_] = histcn([theta_NL],azimuth_edges,'AccumData',dn,'Fun',@sum);
  [dn_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dn,'Fun',@sum);
  [dv_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dv,'Fun',@sum);
  [df_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',df,'Fun',@sum);
  % for iloc = 1:nAzimuth    
  %   dv_tot_all(it,iloc) = sum(MP(it).dv(loc==iloc));
  %   df_tot_all(it,iloc) = sum(MP(it).df(loc==iloc));
  % end
  dn_tot_all_(it,:) = dn_tot_;

  dn_tot_all(it,:,:) = dn_tot;
  dv_tot_all(it,:,:) = dv_tot;
  df_tot_all(it,:,:) = df_tot;
  
end

% Not seperated by E
d3v_tot = sum(PD.d3v('mat'),[2:4]);
d3v_tot_az_bin_ = repmat(d3v_tot/nAzimuth,[1,nAzimuth]);
f_tot_2_ = dn_tot_all_./d3v_tot_az_bin_;

% Seperated by E
d3v_tot_per_E = sum(PD.d3v('mat'),[3:4]);
d3v_tot_az_bin = repmat(d3v_tot_per_E/nAzimuth,[1,1,nAzimuth]);
%f_tot_1 = dn_tot_all./dv_tot_all;
f_tot = dn_tot_all./d3v_tot_az_bin;


iAzim_ = PDist(PD.time,f_tot,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
iAzim_.ancillary.v_shift = vL_shift;
%iAzim_.units = 's^3/m^6';
eval(sprintf('iAzim_%04.0f = iAzim_;',abs(vL_shift)));
%iAzim_ = PDist(PD.time,f_tot_1,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6

%% Calculate some quantitative boundaries
% Divide into top and bottom
isep = nAzimuth/2;


f_tot_use = circshift(iAzim_.elim([1000 Inf]).data,-1,3);
th_tot = squeeze(sum(f_tot_use(:,:,:),2));
th_bot = squeeze(sum(f_tot_use(:,:,1:isep),2));
th_top = squeeze(sum(f_tot_use(:,:,isep+1:nAzimuth),2));
th_bot_cumsum = cumsum(squeeze(sum(f_tot_use(:,:,1:isep),2)),2);
th_top_cumsum = cumsum(squeeze(sum(f_tot_use(:,:,isep+1:nAzimuth),2)),2);

th_bot_cumsum_norm = th_bot_cumsum./repmat(th_bot_cumsum(:,end),1,size(th_bot_cumsum,2));
th_top_cumsum_norm = th_top_cumsum./repmat(th_top_cumsum(:,end),1,size(th_top_cumsum,2));

th_bot_cumsum_norm = smooth2(double(th_bot_cumsum_norm),2,0);
th_top_cumsum_norm = smooth2(double(th_top_cumsum_norm),2,0);

nt = size(f_tot,1);
thresh = 0.2;
th_thresh_bot = zeros(nt,2);
th_thresh_top = zeros(nt,2);
for it = 1:size(f_tot,1)

  i1_bot = find(th_bot_cumsum_norm(it,:)>thresh,1,'first');
  i2_bot = find(th_bot_cumsum_norm(it,:)>(1-thresh),1,'first');
  i1_top = find(th_top_cumsum_norm(it,:)>thresh,1,'first');
  i2_top = find(th_top_cumsum_norm(it,:)>(1-thresh),1,'first');

  th_thresh_bot(it,1) = i1_bot;
  th_thresh_bot(it,2) = i2_bot;
  th_thresh_top(it,1) = i1_top;
  th_thresh_top(it,2) = i2_top;
end

th_thresh_bot = movmean(th_thresh_bot,1,1);
th_thresh_top = movmean(th_thresh_top,1,1);

its = 1:nt;

h = setup_subplots(5,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,log10(th_tot)'); shading(hca,'flat')


hca = h(isub); isub = isub + 1;
pcolor(hca,log10(th_top)'); shading(hca,'flat')
hold(hca,'on')
plot(hca,its,th_thresh_top,'k')
hold(hca,'off')

hca = h(isub); isub = isub + 1;
pcolor(hca,th_top_cumsum_norm'); shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,log10(th_bot)'); shading(hca,'flat')
hold(hca,'on')
plot(hca,its,th_thresh_bot,'k')
hold(hca,'off')

hca = h(isub); isub = isub + 1;
pcolor(hca,th_bot_cumsum_norm'); shading(hca,'flat')

%% Plot
tint_plot = time_xline_ion + 25*[-1 1];
elim = [1500 Inf];
fontsize = 12;
h = irf_plot(2);


if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint_plot),mvaVi?.y.tlim(tint_plot),mvaVi?.z.tlim(tint_plot)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % iAzim_000
  hca = irf_panel('azim iAzim_000');
  specrec = iAzim_0000.elim(elim).tlim(tint_plot).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint_plot,90*[1 1]),'k--')
  irf_plot(hca,irf.ts_scalar(tint_plot,-90*[1 1]),'k--')
  hold(hca,'off')  
  hca.YLabel.String = '\theta_{LN} (deg)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{sprintf('v_L -> v_L-(%g km/s)',iAzim_0000.ancillary.v_shift)},[0.98 0.98],'color','k','fontsize',12)
end
if 1 % iAzim_170
  hca = irf_panel('azim iAzim_170');
  specrec = iAzim_0170.elim(elim).tlim(tint_plot).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel'); 
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint_plot,90*[1 1]),'k--')
  irf_plot(hca,irf.ts_scalar(tint_plot,-90*[1 1]),'k--')
  hold(hca,'off')  
  hca.YLabel.String = '\theta_{LN} (deg)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{sprintf('v_L -> v_L-(%g km/s)',iAzim_0170.ancillary.v_shift)},[0.98 0.98],'color','k','fontsize',12)
end
if 0 % iAzim_999
  hca = irf_panel('azim iAzim_999');
  specrec = iAzim_0999.elim(elim).tlim(tint_plot).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint_plot,90*[1 1]),'k--')
  irf_plot(hca,irf.ts_scalar(tint_plot,-90*[1 1]),'k--')
  hold(hca,'off')  
  hca.YLabel.String = '\theta_{LN} (deg)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{sprintf('v_L -> v_L-(%g km/s)',iAzim_0999.ancillary.v_shift)},[0.98 0.98],'color','k','fontsize',12)
end

h(end).XTickLabelRotation = 0;
colormap('parula')
%colormap(pic_colors('candy4'))
%colormap(flipdim(pic_colors('thermal'),1))
irf_plot_axis_align(h)
irf_zoom(h,'x',tint_plot)

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).CLim = [-30 -26];',1:numel(h))
h(1).Title.String = sprintf('L = [%.2f,%.2f,%.2f]; M = [%.2f,%.2f,%.2f]; N = [%.2f,%.2f,%.2f];',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));

%% Make plot of opening angle based on injection angle

theta = 20;

y = @(x,theta) x*tand(theta);
xx = linspace(0,2,100);
hca = subplot(1,1,1);

plot(hca,xx,y(xx,theta),xx,-y(xx,theta),'k')
hold(hca,'on')
% y = 0 line
plot(hca,xx,0*y(xx,10),'k--')
% quiver of injection directions
x0 = 1; 
y0 = y(x0,theta);
vx = 0.1*sind(theta);
vy = -0.1*cosd(theta);
quiver(hca,x0,y0,vx,vy,0,'k')

y0 = -y(x0,theta);
vx = 0.1*sind(theta);
vy = 0.1*cosd(theta);
quiver(hca,x0,y0,vx,vy,0,'k')

% Add angle info
x1 = 1.5;
x2 = 1.5;
y1 = 0;
y2 = y(x2,theta);
%annotation('doublearrow',[x1 x2],[y1 y2])
text(1.3,0.1,0,sprintf('theta = %g deg.',theta))
%text(1.1,-0.1,0,sprintf('x/y = %.3f/%.3f = %.3f',xx(end),y(xx(end),theta),xx(end)/y(xx(end),theta)))
text(1.1,-0.1,0,sprintf('y/x = %.3f/%.3f = %.3f',y(xx(end),theta),xx(end),y(xx(end),theta)/xx(end)))

hold(hca,'off')

hca.Visible = 'off';
axis(hca,'equal')

%% Plot for paper

tint_plot = [iAzim_0170.time.start iAzim_0170.time.stop];
times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];


time = EpochTT(times_utc(2,:)); % Time of example distribution

h1 = irf_plot(2);

h1(1).Position = [0.5 0.75 0.4 0.15];
h1(2).Position = [0.5 0.15 0.4 0.6];

h2 = subplot(1,2,1);
h2.Position = [0.1 0.15 0.5 0.5];
h2.Position = [0.1181    0.2759    0.2267    0.5810];


if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint_plot),mvaVi?.y.tlim(tint_plot),mvaVi?.z.tlim(tint_plot)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % iAzim_170
  hca = irf_panel('azim iAzim_170');
  specrec = iAzim_0170.elim(elim).tlim(tint_plot).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel'); 
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint_plot,90*[1 1]),'k--')
  irf_plot(hca,irf.ts_scalar(tint_plot,-90*[1 1]),'k--')
  hold(hca,'off')  
  hca.YLabel.String = '\theta_{LN} (deg)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{sprintf('v_L -> v_L-(%g km/s)',iAzim_0170.ancillary.v_shift)},[0.98 0.98],'color','k','fontsize',12)
  hca.XGrid = 'on'; hca.YGrid = 'on';
  hca.CLim = [-30 -26];
end

irf_plot_axis_align(h1)
h1(1).YLim(1) = [-799];



elim = [000 Inf];
%elim = [200 Inf];
%c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)
%c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)
pdist_all = pdist_cleaned;

  pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
  tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
  elow = max(tsElow.tlim(tint_dist).data);  
  %elow = 2000;
  pdist = pdist.elim([elow Inf]);
  %pdist.data(:,:,:,1) = 0;
  %pdist.data(:,:,:,end) = 0;

  % Calculate eigenvectors and values
  pmoms_tmp = mms.psd_moments(pdist,scPot3);
  mvaP = irf.ts_tensor_xyz(pmoms_tmp.P_psd.time,lmn*squeeze(pmoms_tmp.P_psd.data)*lmn');
  [tsEig_val_2d_tmp, tsEig_v1_2d_tmp, tsEig_v2_2d_tmp] = mvaP.eig([1 2]);

  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))
  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;

  isub = 1;
  if 1 % f(L,N)
    hca = h2(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl],'vint',vint_M);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 3500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.XTick = -4000:2000:4000;
    hca.YTick = -4000:2000:4000;
    hca.XMinorGrid = 'on';
    hca.YMinorGrid = 'on';
    %hca.XMinorTick = -4000:2000:4000;
    %hca.YMinorTick = -3000:1000:3000;
    if 1 % add circular grid
      %%
      hold(hca,'on')      
      ang = linspace(0,360,361);
      r = 2200;
      color_grid = ([0.9 0.9 0.9]-0.3)*0;
      x0 = vL_shift;
      y0 = 0;
      %plot(hca,x0 + r*cosd(ang), y0 + r*sind(ang),'color',color_grid)
      plot(hca,x0 + 1500*cosd(ang), y0 + 1500*sind(ang),'color',color_grid)
      %plot(hca,x0 + 2000*cosd(ang), y0 + 2000*sind(ang),'color',color_grid)
      plot(hca,x0,0,'kx',x0,0,'ko','linewidth',2)
      angs = -180:30:180;
      for ia = 1:numel(angs)
        xx = x0 + [0 r*cosd(angs(ia))];
        yy = y0 + [0 r*sind(angs(ia))];
        if abs(angs(ia)) == 90
          plot(hca,xx,yy,'color',color_grid,'linestyle','--','linewidth',2)
        else
          plot(hca,xx,yy,'color',color_grid,'linestyle','-')

        end
        text_margin = 1.3;
        if angs(ia) == 180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%4.f^o',angs(ia))}, ...
            'color','k','VerticalAlignment','bottom','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        elseif angs(ia) == -180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%g^o',angs(ia))}, ...
            'color','k','VerticalAlignment','top','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        else
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,sprintf('%4.f^o',angs(ia)), ...
            'color','k','VerticalAlignment','middle','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        end

      end
      
      hold(hca,'off')
    end
  end


h1(end).XTickLabelRotation = 0;
%colormap('parula')
cmap = irf_colormap('magma');
colormap(cmap)
%colormap(flipdim(irf_colormap('Spectral'),1))
%irf_colormap('vik')
%colormap(pic_colors('candy4'))
%colormap(flipdim(pic_colors('thermal'),1))

irf_zoom(h1,'x',tint_plot)

irf_legend(h2,'(a)',[-0.05 1.05],'color','k','fontsize',16)
irf_legend(h1(1),'(b)',[-0.1 1.01],'color','k','fontsize',16)
irf_legend(h1(2),'(c)',[-0.1 1.0],'color','k','fontsize',16)

%h(1).Title.String = sprintf('L = [%.2f,%.2f,%.2f]; M = [%.2f,%.2f,%.2f]; N = [%.2f,%.2f,%.2f];',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));


%%

R = linspace(1,10,100);
  theta = linspace(-180,180,360);
  Z = linspace(0,10,360)'*linspace(0,10,100);
  figure
  polarPcolor(R,theta,Z*nan,'Ncircles',3)
  hold off