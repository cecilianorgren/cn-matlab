%% Get list of magnetosheath time intervals
localuser = 'cecilia';
localuser = 'cno062';
table_ubmshc = readtable(['/Users/' localuser '/MATLAB/cn-matlab/projects/plasma_heating/data/mms_unbiased_magnetosheath_campaign.txt']);
region_prob = load(['/Users/' localuser '/Data/MMS/DB_Lalti/Proba_full.mat']);

%% Load data, unbiased magnetosheath campaign
%mms.db_init('local_file_db','/Volumes/mms');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS/']);
db_info = datastore('mms_db');

nIntervals = size(table_ubmshc,1);

files = struct([]);

for iInterval = 1:nIntervals
  t1 = table_ubmshc.time{iInterval}; t1 = [t1(1:10), 'T', t1(12:19), '.0Z'];
  t2 = table_ubmshc.xEnd{iInterval}; t2 = [t2(1:10), 'T', t2(12:19), '.0Z'];
  tint = EpochTT([t1; t2]);
  
  files_tmp = mms.db_list_files('mms1_fgm_brst_l2',tint);
  files = cat(1,files,files_tmp');
end

%% Get local file list
tint_all = irf.tint('2015-01-01T00:00:00.00Z/2024-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

%% Go through burst intervals to identify current sheets

nFiles = size(files,1);
for iFile = 1%nFiles-3%:nFiles
  tint = [files(iFile).start files(iFile).stop];
  
  % Load data  
  c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
  c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
  [Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');


  %% Find current sheets
  % Box average current to given window
  dT = 10;
  box_timeline = tint(1):dT:tint(2);
  Jcurl_box = Jcurl.resample(box_timeline);

  deltaJ = Jcurl + -1*Jcurl_box.resample(Jcurl);
  
  % Criteria
  Threshold = 10e-9;
  MinPeakProminence = 20e-9;
  [PKS,LOCS,W] = findpeaks(Jcurl.abs.data,'MinPeakProminence',MinPeakProminence);

  % Plot data
  
  colors = pic_colors('matlab');
  h = irf_plot(4);
  
  hca = irf_panel('B');
  set(hca,'colororder',mms_colors('xyza'))
  irf_plot(hca,{gseB.x,gseB.y,gseB.z,gseB.abs},'comp')
  
  hca = irf_panel('Jcurl');
  set(hca,'colororder',mms_colors('xyza'))
  irf_plot(hca,{Jcurl.x,Jcurl.y,Jcurl.z,Jcurl.abs},'comp')

  hca = irf_panel('Jabs');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{Jcurl.abs,Jcurl_box.abs},'comp')
  
  hca = irf_panel('delta J');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{deltaJ.abs},'comp')

  % Plot peaks
  time_pks = Jcurl.time(LOCS);  
  %tint_mark = [time_pks.epochUnix-0.5*W, time_pks.epochUnix+0.5*W];
  tint_mark = time_pks.epochUnix;
  for iTint = 1:numel(LOCS)
    irf_pl_mark(h,tint_mark(iTint,:),'k')
  end

end

%% Get the Gaussian Mixture Model for a FPI distribution
%c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
%c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
units = irf_units;

Threshold = 10e-9;
MinPeakProminence = 20e-9;
[PKS,LOCS,W] = findpeaks(Jcurl.abs.data,'MinPeakProminence',MinPeakProminence);
iLargestPeak = LOCS(find(max(PKS)));

pdist = ePDist1(iLargestPeak);

% To use a Gaussian Mixture Model, which is a non-weighted method, the
% macroparticles should be many enough so that they have a very similar
% weight.
nParticlesTarget = 32*32*16*50;
mParticles = pdist.macroparticles('ntot',nParticlesTarget,'skipzero',1);
% hist(macroParticles.dn) 

% Data which to fit the GMM to
R = [mParticles.vx, mParticles.vy, mParticles.vz];

% Define a grid on which to plot the results
Emax = 10e3; % eV
vmax = sqrt(Emax*units.eV*2/units.me)*1e-3; % km/s
xvec = linspace(-vmax,vmax,101);
yvec = linspace(-vmax,vmax,102);
zvec = linspace(-vmax,vmax,103);
[X,Y,Z] = ndgrid(xvec,yvec,zvec);
XYZ = [X(:) Y(:) Z(:)];

% Fit model
nComp = 3;
gm = fitgmdist(R,nComp);

% The pdf function only works 1D data, so to be able to plot N-dimensional
% data, arrayfun is used to apply it to each 'column' (row?) of the
% higher-dimensional array
gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm,[x0 y0 z0]),x,y,z);
F = gmPDF(X,Y,Z);  

F1 = gm.ComponentProportion(1)*mvnpdf(XYZ, gm.mu(1,:), gm.Sigma(:,:,1)); F1 = reshape(F1,size(X));
F2 = gm.ComponentProportion(2)*mvnpdf(XYZ, gm.mu(2,:), gm.Sigma(:,:,2)); F2 = reshape(F2,size(X));
F3 = gm.ComponentProportion(3)*mvnpdf(XYZ, gm.mu(3,:), gm.Sigma(:,:,3)); F3 = reshape(F3,size(X));

% Get index of which cluster each data ount belong to, for plotting
clusterR = cluster(gm,R);

%
if 0
  %%
  tic;
  eva = evalclusters(R,'gmdistribution','silhouette','KList',1:3);
  toc
  
end

% Normalize to f based on clustering
n = sum(mParticles.dn);
n1 = sum(mParticles.dn)*gm.ComponentProportion(1);
n2 = sum(mParticles.dn)*gm.ComponentProportion(2);
n3 = sum(mParticles.dn)*gm.ComponentProportion(3);

if 0
f1 = sum(mParticles.df)*gm.ComponentProportion(1);
f2 = sum(mParticles.df)*gm.ComponentProportion(2);
f3 = sum(mParticles.df)*gm.ComponentProportion(3);

f1_b = sum(mParticles.df(clusterR==1));
f2_b = sum(mParticles.df(clusterR==2));
f3_b = sum(mParticles.df(clusterR==3));

f1_c = mean(mParticles.df(clusterR==1));
f2_c = mean(mParticles.df(clusterR==2));
f3_c = mean(mParticles.df(clusterR==3));
end
%n = sum(mParticles.dn);
%n1 = sum(mParticles.dn)*gm.ComponentProportion(1);
%n2 = sum(mParticles.dn)*gm.ComponentProportion(2);
%F3 = F3*n3;

%%
nrows = 4;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Reduced 1D, x
  hca = h(isub); isub = isub + 1;  
  vdf = pdist.reduce('1D',[1 0 0]);
  plot(hca,vdf.depend{1}*1e-3,vdf.data);
  if 0
    hold(hca,'on')
    dvy = (yvec(2) - yvec(1))*1e3; % m/s
    dvz = (zvec(2) - zvec(1))*1e3; % m/s
    pl = squeeze(sum(sum(F,3),2))*dvy*dvz;
    plot(hca,xvec,pl)
    hold(hca,'off')
  end
end
if 1 % Reduced 1D, x
  hca = h(isub); isub = isub + 1;  
  dvy = (yvec(2) - yvec(1))*1e3; % m/s
  dvz = (zvec(2) - zvec(1))*1e3; % m/s
  pl = squeeze(sum(sum(F,3),2))*dvy*dvz;

  plot(hca,xvec*1e-3,pl)  
end

if 1 % Reduced 2D, xy
  hca = h(isub); isub = isub + 1;  
  vdf = pdist.reduce('2D',[1 0 0],[0 1 0]);
  vdf.plot_plane(hca);
end
if 1 % Reduced 2D, xx
  hca = h(isub); isub = isub + 1;  
  vdf = pdist.reduce('2D',[1 0 0],[0 0 1]);
  vdf.plot_plane(hca);
end
if 1 % Reduced 2D, fitted distribution 
  hca = h(isub); isub = isub + 1;    
  %%
  dv = (zvec(2) - zvec(1))*1e3; % m/s
  pl = squeeze(sum(F,3))*dv;
  pcolor(hca,xvec*1e-3,yvec*1e-3,log10(pl)')
  shading(hca,'flat')
  hcb = colorbar(hca);
end
if 1 % Iso surface plot of data
  hca = h(isub); isub = isub + 1;
  nSmooth = 3;
  %hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);
  hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill');
  
  c_eval('hs.Patch(?).FaceAlpha = 0.5;',1:numel(hs.Patch))
end
if 1 % Scatter plot of macroparticles with isosurfaces
  hca = h(isub); isub = isub + 1;
  
  nScatterPoints = 1e3;
  iPlot = ceil(linspace(2,numel(mParticles.vx)-1,nScatterPoints));
  scatter3(hca,R(iPlot,1),R(iPlot,2),R(iPlot,3),20,clusterR(iPlot),'marker','.') % Scatter plot with points of size 10
  
  hold(hca,'on')
    
  
  % Flev = iso_values(isurf);
  valIso = 2*1e-3;
  s = isosurface(X,Y,Z,F);
  s1 = isosurface(X,Y,Z,F1);
  s2 = isosurface(X,Y,Z,F2);
  s3 = isosurface(X,Y,Z,F3);
  
  %p = patch(hca,'Faces',s.faces,'Vertices',s.vertices,'FaceAlpha',0.5);
  p1 = patch(hca,'Faces',s1.faces,'Vertices',s1.vertices,'FaceAlpha',0.2,'FaceColor','r','EdgeColor','none');
  p2 = patch(hca,'Faces',s2.faces,'Vertices',s2.vertices,'FaceAlpha',0.2,'FaceColor','b','EdgeColor','none');
  p3 = patch(hca,'Faces',s3.faces,'Vertices',s3.vertices,'FaceAlpha',0.2,'FaceColor','g','EdgeColor','none');
  
  lighting(hca,'gouraud')
  camlight(hca);
  
  hold(hca,'off')
end

