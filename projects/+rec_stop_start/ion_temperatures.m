%% Load data
% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic_dist = 1;
ic = 1:4;

localuser = datastore('local','user');
localuser = 'cecilia';
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
%mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
db_info = datastore('mms_db'); 

tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% Time from time interval
tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');
tint_action = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');
tint = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');

% Time from file name
fileId = '20170725220853';

iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));

% Load data
tic
disp('Loading data.')
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',1:4);
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
disp('Loading VDFs.')
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic_dist) % missing some ancillary data
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic_dist) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?_counts.data./iPDistErr?.data).^2;',ic_dist)
% counts = (data/err)^2 --> 1count = sqrt(1)*err;
c_eval('iPDist?_onecount = iPDistErr?;',ic_dist)


toc
disp('Done.')

%% Prepare data
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)
c_eval('Ti?perp = 0.5*(facTi?.yy+facTi?.zz);',ic)
c_eval('Ti?par = facTi?.xx;',ic)

c_eval('iPDist?_cleaned = iPDist?.movmean(4,''RemoveOneCounts'',iPDist?_counts);',ic_dist)

elim = [700 Inf];
vint = [-Inf Inf];
c_eval('if1Dx? = iPDist?.elim(elim).reduce(''1D'',[1 0 0],''vint'',vint);',ic_dist)
c_eval('if1Dx?_cleaned = iPDist?_cleaned.reduce(''1D'',[1 0 0],''vint'',vint);',ic_dist)

%% Estimate temperatures

[h1,h2] = initialize_combined_plot('leftright',5,2,2,0.4,'vertical');

hca = irf_panel('B');
hca.ColorOrder = mms_colors('xyz1');
irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z,gseB1.abs},'comp')
hca.YLabel.String = 'B_{GSE} (nT/m)';

hca = irf_panel('n');
irf_plot(hca,{ne1},'comp')
hca.YLabel.String = 'n_e (cc)';

hca = irf_panel('Ti');
irf_plot(hca,{(gseTi1.xx+gseTi1.yy+gseTi1.zz)/3},'comp')
hca.YLabel.String = 'T_i (eV)';

hca = irf_panel('deflux');
irf_spectrogram(hca,iPDist1.omni.deflux.specrec)
colormap(hca,pic_colors('pasteljet'))
hca.YScale = 'log';

if 0 
  hca = irf_panel('deflux cleaned');
  irf_spectrogram(hca,iPDist1_cleaned.omni.deflux.specrec)
  colormap(hca,pic_colors('pasteljet'))
  hca.YScale = 'log';
end

hca = irf_panel('fi(vx)');
irf_spectrogram(hca,if1Dx1.specrec('velocity'))
colormap(hca,pic_colors('pasteljet'))

if 0
  hca = irf_panel('fi(vx) cleaned');
  irf_spectrogram(hca,if1Dx1_cleaned.specrec('velocity'))
  colormap(hca,pic_colors('pasteljet'))
end

h1(end).XTickLabelRotation = 0;
%irf_plot_axis_align(h1)
c_eval('h1(?).Position(3) = h1(end).Position(3);',1:numel(h1))
irf_zoom(h1,'x',gseB1.time)



times_dist = EpochTT(['2017-07-25T22:10:03.00Z';...
                      '2017-07-25T22:10:07.50Z']);

%times_dist = EpochTT(['2017-07-25T22:09:55.00Z';...
%                      '2017-07-25T22:10:07.50Z']);
                    
times_dist_exact = EpochTT([]);
clear times_dist_exact_eu hmarkb hmark_tmp times_dist_exact_eu_edges

elim = [500 Inf];

ic = 1;
R = [1 0 0; 0 1 0; 0 0 1];
ip = 0;
ndists = 6;
%fmax = @() n/(2*pi)^0.5

e = 1.6022e-19;
kB = 1.38e-23;
me = 9.1094e-31;
m = 9.1094e-31;
mp = 1.6726e-27;
w = @(T) sqrt(2*units.eV*T/mp);
fmax = @(v,T,n,vd) n/((pi)^(1/2)*w(T))*exp(-(v-vd).^2./w(T)./w(T));

% 2D
for id = 1:times_dist.length
  ip = ip + 1;
  hca = h2(ip);
  dt = 0.05;
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  v1 = R(1,:);
  v2 = R(2,:);
  v3 = R(3,:);
  x_str = 'v_x (km/s)';
  y_str = 'v_y (km/s)';
  z_str = 'v_z (km/s)';
  % bugcheck, slightly rotate system
   % rotang = 30;
   % v2n = v2*cosd(rotang) + v3*sind(rotang);
   % v3n = -v2*sind(rotang) + v3*cosd(rotang);
   % v2 = v2n; v3 = v3n;
  vlim = 2300;
  vg = -vlim:100:vlim;
  if 1 % x,z
    
    %c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+6*0.15*0.5*[-1 1]+0);',ic)
    c_eval('ddist = iPDist?.elim(elim).movmean(ndists).tlim(times_dist(id)+1*0.15*0.5*[-1 1]);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = mean(ddist.time.epochUnix);
    times_dist_exact_eu_edges(id,1) = ddist.time(1).epochUnix-0.03*0.5;
    times_dist_exact_eu_edges(id,2) = ddist.time(end).epochUnix+0.03*0.5;
    f2d = ddist.elim([00 Inf]).reduce('2D',v1,v3);

    f2d.plot_plane(hca,'smooth',0);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = x_str;
    hca.YLabel.String = z_str;
    
    hmark_tmp = irf_pl_mark(h1,ddist.time+ndists*0.5*0.03*[-1 1],'k','facealpha',0.1);
    drawnow
    %c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

  end
  if 0 % x,y
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+0.15*0.5*[-1 1]);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vg',vg);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  colormap(hca,pic_colors('pasteljet'))
end

% 1D
ip = ip + 1;  
for id = 1:times_dist.length  
%  ip = ip + 1;
  hca = h2(ip);
  hca.NextPlot = "add";
  hca.ColorOrder = pic_colors('matlab');
  dt = 0.05;
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  v1 = R(1,:);
  v2 = R(2,:);
  v3 = R(3,:);
  x_str = 'v_x (km/s)';
  y_str = 'v_y (km/s)';
  z_str = 'v_z (km/s)';

  vlim = 2300;
  vg = -vlim:100:vlim;
  %elim = [200 Inf];
  if 1 % x
    
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+6*0.15*0.5*[-1 1]+0);',ic)
    %c_eval('ddist = iPDist?.elim(elim).movmean(ndists).tlim(times_dist(id)+1*0.15*0.5*[-1 1]);',ic)
    c_eval('ddist_1c = iPDist?_onecount.elim(elim).tlim(times_dist(id)+6*0.15*0.5*[-1 1]+0);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = mean(ddist.time.epochUnix);
    times_dist_exact_eu_edges(id,1) = ddist.time(1).epochUnix-0.03*0.5;
    times_dist_exact_eu_edges(id,2) = ddist.time(end).epochUnix+0.03*0.5;
    f1d = ddist.elim(elim).reduce('1D',v1);
    f1d_counts = ddist.elim(elim).reduce('1D',v1,'counts',1);
    f1d_onecount = ddist_1c.elim(elim).reduce('1D',v1,'counts',1);

   % plot(hca,f1d.depend{1}(1,:),0*mean(f1d.data,1),...
   %   f1d_counts.depend{1}(1,:),mean(f1d_counts.data,1),'--',...
   %   f1d_onecount.depend{1}(1,:),mean(f1d_onecount.data,1),'-.')

   if id == 1
    plot(hca,f1d.depend{1}(1,:),mean(f1d.data,1),...
      f1d_counts.depend{1}(1,:),mean(f1d_counts.data,1),'--','linewidth',2)
   elseif id == 2
    plot(hca,f1d.depend{1}(1,:),mean(f1d.data,1)*1,...
      f1d_counts.depend{1}(1,:),mean(f1d_counts.data,1),'--','linewidth',2)
   end

    if 0
      %%
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,10000,0.16*1e6,800e3),'k--')
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,6000,0.09*1e6,500e3),'k--')
      %plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,6000,0.08*1e6,500e3)...
      %                           +fmax(f1d.depend{1}(1,:)*1e3,600,0.01*1e6,500e3),'r-')

    elseif 0
      %%
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,11000,0.2*1e6,700e3),'k:')
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,6000,0.1*1e6,500e3),'k:')
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,6000,0.08*1e6,500e3)...
                                 +fmax(f1d.depend{1}(1,:)*1e3,600,0.01*1e6,500e3),'r:')
    elseif 0
      %%
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,8000,0.15*1e6,500e3),'k--')
      plot(hca,f1d.depend{1}(1,:),fmax(f1d.depend{1}(1,:)*1e3,5000,0.16*1e6,500e3),'k--')
    end
    hca.YScale = 'lin';
    hca.Layer = 'top';
    axis(hca,'square')
    hca.XLabel.String = x_str;
    %hca.YLabel.String = z_str;
    hca.YLabel.String = ['f_i' f1d.units];
    
    hmark_tmp = irf_pl_mark(h1,ddist.time,'k','facealpha',0.1);
    drawnow
    %c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

  end
 
  colormap(hca,pic_colors('pasteljet'))
end

%hlinks = linkprop(h2,{'CLim','XLim','YLim'});
%h2(1).CLim = [-10 -7.];
