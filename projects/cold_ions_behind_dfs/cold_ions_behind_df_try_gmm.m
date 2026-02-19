%% Load DB table
localuser = 'cecilianorgren';
%file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2022.nc'];

file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.nc'];
%data = load(file);
file_csv = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.csv'];
tlim = readtable(file_csv);
tstart = tlim.Var1; tstart.TimeZone = "UTCLeapSeconds";
tstop = tlim.Var2;  tstop.TimeZone = "UTCLeapSeconds";

tstart_ttns = convertTo(tstart, 'tt2000');
tstop_ttns = convertTo(tstop, 'tt2000');
t_duration_s = double((tstop_ttns - tstart_ttns))*1e-9;

%ncdisp(file)

info = ncinfo(file);
vars = {info.Variables.Name};
nvars = numel(vars);

clear db
for ivar = 1:nvars
  db.(vars{ivar}) = ncread(file,vars{ivar});
end

db_table_ff = struct2table(db);
% Add stop time of slow
db_table_ff = addvars(db_table_ff,tstart_ttns,tstop_ttns,t_duration_s,'After','time');


t0_ff = EpochTT('2017-05-05T19:41:44.790324');
% time: microseconds since 2017-05-05T19:41:44.790324
time_ff_ttns = t0_ff.ttns + db_table_ff.time*1e3;
db_table_ff.time = time_ff_ttns; % rewrite time in ttns


t0_df = EpochTT('2017-05-19T03:06:44.458185978');
% t_df: nanoseconds since  2017-05-19T03:06:44.458185978
time_df_ttns = double(t0_df.ttns + int64(db_table_ff.t_df));
time_df_ttns(~db_table_ff.is_df) = NaN;
db_table_ff.t_df = time_df_ttns; % rewrite time in ttns

%% Prepare table
if 0
clear table_ci
gm_before = {};
gm_after = {};
T_scalar_before = {};
T_scalar_after = {};
T_before = {};
T_after = {};
T_min_before = {};
T_min_after = {};
id_df = [];
t_ff_start = string([]);
t_ff_stop = string([]);
t_df = string([]);
table_ci = table(id_df,t_ff_start,t_ff_stop,t_df,gm_before,gm_after,T_before,T_after,T_scalar_before,T_scalar_after,T_min_before,T_min_after);
end
%% Loop events and plot brief overview
units = irf_units;
ic = 1;


db_table_df = db_table_ff(db_table_ff.is_df==1,:);
nDF = numel(db_table_df.time);


iDFs = [9 128];
iDFs = 217:nDF;
%iDFs = 87;
iDFs = 17;
doPrint = 0;
doPlot = 1;
for iDF = iDFs%87%iDFs(1)
  try
  disp(iDF)
  % Define time
  t0 = EpochTT(db_table_df.time(iDF));
  T =  db_table_df.t_duration_s(iDF);  
  t_ff_start = t0;
  t_ff_stop = t_ff_start + T;
  tint = [t_ff_start t_ff_stop] + 5*[-1 1];
  tDF = EpochTT(int64(db_table_df.t_df(iDF)));

  % Load data
  c_eval('gseB = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  if isempty(gseB); disp('Skipping.'); continue; end
  c_eval('gseVi = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
  if isempty(gseVi); disp('Skipping.'); continue; end
  c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
  c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
  %c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
  %iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;

  %% Prep data
  %PD_use = iPDist_counts; % iPDist
  PD_use = iPDist; % iPDist
  tsElow = PD_use.find_noise_energy_limit(5).movmean(15);
  emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit
  PD = PD_use.mask({emask_mat});
  %nMovMean = 5; % Number of distributions for the moving average
  %PD = PD.movmean(3);

  %% Do the Gaussian Mixture Model
  nMP = 200000; % Number of macroparticles
  nGroupsMax = 4;  % number of classes/groups for kmeans and gmm
  tShift = 5; % apply to two times, before and after DF

  tPreDF  = tDF + -tShift;
  tPostDF = tDF + +tShift;
  times = [tPreDF tPostDF];

  % Initial guess
  initial_guess = [0 -1000 -1000;...
                   0 -1000 +1000;...
                   0 +1000 -1000;...
                   0 +1000 +1000;...
                   0  0000  0000;...
                   0  0000  0000]*0.1;


  nt = times.length;
  %clear gm ntot T T_scalar;
  
%
  for it = 1:nt
    time = times(it);
    pdist = PD.tlim(time+0.5*0.150*[-1 1]);
    %nMP = nansum(round(pdist.data(~isnan(pdist.data)))); % total number of counts
    
    scpot = mean(scPot.tlim(time + 0.5*0.15*[-1 1]).data,1);
    scpot = irf.ts_scalar(time,scpot);
    B = mean(gseB.tlim(time + 0.5*0.15*[-1 1]).data,1); B = B/norm(B);
    
    for nGroups = 1:nGroupsMax
      MP = pdist.macroparticles('ntot',nMP,'skipzero',1,'scpot',scpot);
      nMP = numel(MP.dv);
      MP.dn = MP.df.*MP.dv;
      V_dbcs = [MP.vx, MP.vy, MP.vz]; 
      
      % Need to rotate these into the specified coordinate system  
      % The rotation is not done exactly right, due to not taking into account
      % differences between dbcs and gse. To do in the future
      %V = V_dbcs*lmn';
      V = V_dbcs;
      MP.vx = V(:,1);
      MP.vy = V(:,2);
      MP.vz = V(:,3);
      %MP
          
      S.mu = initial_guess(1:nGroups,:);
      S.Sigma = repmat([500 10 10; 10 500 10; 10 10 500],[1 1 nGroups]);  
      S.ComponentProportion = repmat(1,[1,nGroups]);

      options = statset('Display','final','MaxIter',1500,'TolFun',1e-5);


      X = [MP.vx, MP.vy, MP.vz];
      %X(sqrt(sum(X.^2,2))>1000,:) = [];
      %gm{it,nGroups} = fitgmdist(X,nGroups,'Start',S,'SharedCovariance',false,'Options',options);
      gm{it,nGroups} = fitgmdist(X,nGroups,'Start','randSample','SharedCovariance',false,'Options',options,'RegularizationValue',0.1);

       
      %X = [MP.vx, MP.vy, MP.vz];   
      %S.mu = gm{it,nGroups}.mu;
      %S.Sigma = gm{it,nGroups}.Sigma;
      %S.ComponentProportion = gm{it,nGroups}.ComponentProportion;
      %gm{it,nGroups} = fitgmdist(X,nGroups,'Start','randSample','SharedCovariance',false,'Options',options,'RegularizationValue',0.1);
      
      ntot(it) = sum(MP.df.*MP.dv); % cc;

      % Extract som physical quantities
      vt2 = gm{it,nGroups}.Sigma*1e6; % (km/s)^2 -> (m/s)^2
      cov_T_mat = units.mp*vt2/2/units.eV;
      T_mat{it,nGroups} = cov_T_mat;
      for iGroup = 1:nGroups % Lopp through gaussion components
        T_scalar{it,nGroups}(iGroup) = trace(T_mat{it,nGroups}(:,:,iGroup))/3;
      end
     
    end    
  end

  % Diagnose cold-ion quantities
  % Find coldest temperature for each it, K
  T_min = cellfun(@(x)min(x),T_scalar);
  T_max = cellfun(@(x)max(x),T_scalar);

  % Add results to table
  if 0
  table_ci.id_df(iDF) = iDF;
  table_ci.t_ff_start(iDF) = t_ff_start.utc;
  table_ci.t_ff_stop(iDF) = t_ff_stop.utc;
  table_ci.t_df(iDF) = tDF.utc;
  table_ci.gm_before(iDF) = {gm(1,:)};
  table_ci.gm_after(iDF) ={gm(2,:)};
  table_ci.T_before(iDF) = {T_mat(1,:)};
  table_ci.T_after(iDF) = {T_mat(2,:)};
  table_ci.T_scalar_before(iDF) = {T_scalar(1,:)};
  table_ci.T_scalar_after(iDF) = {T_scalar(2,:)};
  table_ci.T_min_before(iDF) = {T_min(1,:)};
  table_ci.T_min_after(iDF) = {T_min(2,:)};
  table_ci.T_max_before(iDF) = {T_max(1,:)};
  table_ci.T_max_after(iDF) = {T_max(2,:)};
  end


  % Plot results  
  if doPlot
    fontsize = 13;
    vlim = 2500;
    nMC = 500;
  
    xvec = linspace(-vlim,vlim,101); dvx = xvec(2)-xvec(1);
    yvec = linspace(-vlim,vlim,102); dvy = yvec(2)-yvec(1);
    zvec = linspace(-vlim,vlim,103); dvz = zvec(2)-zvec(1);
    [X,Y,Z] = ndgrid(xvec,yvec,zvec);
    XYZ = [X(:) Y(:) Z(:)];
  
    [h1,h2] = initialize_combined_plot('topbottom',3,2,5,0.3,'vertical');
    
    hca = irf_panel('B');
    hca.ColorOrder = mms_colors('xyza');
    irf_plot(hca,{gseB.x,gseB.y,gseB.z},'comp')
    hca.YLabel.String = 'B (nT)';
  
    hca = irf_panel('Vi');
    hca.ColorOrder = mms_colors('xyza');
    irf_plot(hca,{gseVi.x,gseVi.y,gseVi.z},'comp')
    hca.YLabel.String = 'v_i (nT)';
  
    hca = irf_panel('ion deflux omni');
    irf_spectrogram(hca,PD.deflux.omni.specrec,'log')
    hca.YScale = "log";
    hca.Color = [0 0 0] + 0.95;
    
    hmark = irf_pl_mark(h1,tDF,'k','linewidth',1);
    %userdata = get(gcf,'userdata');
    text(h1(1),hmark(1).XData(1),h1(1).YLim(2)*0.99,'DF','VerticalAlignment','bottom', "HorizontalAlignment","center",'FontSize',fontsize)
  
    irf_plot_axis_align
    irf_zoom(h1,'x',[PD.time.start PD.time.stop])
    c_eval('h1(?).YLabel.Interpreter = ''tex'';',1:numel(h1))
    h1(end).XTickLabelRotation = 0;
    c_eval('h1(?).FontSize = fontsize;',1:numel(h1))
    colormap([flipdim(irf_colormap('Spectral'),1)])
  
    for it = 1:nt
      hmark_tmp = irf_pl_mark(h1,times(it),'r');
      text(h1(1),hmark_tmp(1).XData(1),h1(1).YLim(2)*1.5,sprintf('VDF%g',it),'VerticalAlignment','bottom', "HorizontalAlignment","center",'FontSize',fontsize)
    end
    
    isub = 1;
    % 2D vdfs  
    for it = 1:nt
      hca = h2(isub); isub = isub + 1;
      time = times(it);
      pdist = PD.tlim(time+0.5*0.150*[-1 1]);
      vdf = pdist.reduce('2D',[0.99 0 0],[0 0 0.99]);
      vdf.plot_plane(hca)
      hca.Title.String = sprintf('VDF%g',it);
      hca.XLabel.String = 'v_x (km/s)';
      hca.YLabel.String = 'v_z (km/s)';
      axis(hca,'square')
      hca.XLim = vlim*[-1 1];
      hca.YLim = vlim*[-1 1];
      hold(hca,'on')
      quiver(hca,-B(1)*vlim,-B(2)*vlim,B(1)*2*vlim,B(2)*2*vlim,0,'k')
      hold(hca,'off')
    end
    if 1 % Plot fits for each K
      for it = 1:nt
        hca = h2(isub); isub = isub + 1;
        time = times(it);
        pdist = PD.tlim(time+0.5*0.150*[-1 1]);
        vdf = pdist.reduce('1D',[0 0 1],'nMC',nMC);
        plot(hca,vdf.depend{1},vdf.data,'linewidth',1.5)
        hold(hca,'on')
        for K = 1:nGroupsMax
          mu = gm{it,K}.mu;
          Sigma = gm{it,K}.Sigma;
          gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm{it,K},[x0 y0 z0]),x,y,z);
          
          Ftot = X*0;
          for iComp = 1:K
            Ftmp = gm{it,K}.ComponentProportion(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
            F{iComp} = reshape(Ftmp,size(X));  
            Ftot = Ftot + F{iComp};
          end
          Fplot = squeeze(sum(Ftot,[1 2]))*ntot(it)*dvx*dvy*1e3;     
          plot(hca,zvec,Fplot)
        end
        hca.XLabel.String = 'v_z (km/s)';
        hca.YLabel.String = 'f(v_z) (s/m^4)';
        irf_legend(hca,["Observed" "K="+(1:nGroupsMax)]',[0.98 0.98],'fontsize',fontsize)
      end
      hold(hca,'off')
    end
  
    if 1 % Plot component fits for one K
      K = 2;
      for it = 1:nt
        hca = h2(isub); isub = isub + 1;
        time = times(it);
        pdist = PD.tlim(time+0.5*0.150*[-1 1]);
        vdf = pdist.reduce('1D',[0 0 1],'nMC',nMC);
        plot(hca,vdf.depend{1},vdf.data,'linewidth',1.5)
        hold(hca,'on')
        mu = gm{it,K}.mu;
        Sigma = gm{it,K}.Sigma;
        gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm{it,K},[x0 y0 z0]),x,y,z);
        
        Ftot = X*0;
        for iComp = 1:K
          Ftmp = gm{it,K}.ComponentProportion(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
          F{iComp} = reshape(Ftmp,size(X));  
          Ftot = Ftot + F{iComp};
          Fplot = squeeze(sum(F{iComp},[1 2]))*ntot(it)*dvx*dvy*1e3;  
          plot(hca,zvec,Fplot)  
        end
        Fplot = squeeze(sum(Ftot,[1 2]))*ntot(it)*dvx*dvy*1e3;     
        plot(hca,zvec,Fplot)
        hca.XLabel.String = 'v_z (km/s)';
        hca.YLabel.String = 'f(v_z) (s/m^4)';
        irf_legend(hca,["Observed" "iComp = "+(1:K) "All comp."]',[0.98 0.98],'fontsize',fontsize)
        irf_legend(hca,["Observed","T = "+round(T_scalar{it,K})+" eV", "..."]',[0.02 0.98],'fontsize',fontsize)
      end
      hold(hca,'off')
    end
  
    if 1 % Plot component fits for one K
      K = 3;
      for it = 1:nt
        hca = h2(isub); isub = isub + 1;
        time = times(it);
        pdist = PD.tlim(time+0.5*0.150*[-1 1]);
        vdf = pdist.reduce('1D',[0 0 1],'nMC',nMC);
        plot(hca,vdf.depend{1},vdf.data,'linewidth',1.5)
        hold(hca,'on')
        mu = gm{it,K}.mu;
        Sigma = gm{it,K}.Sigma;
        gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm{it,K},[x0 y0 z0]),x,y,z);
        
        Ftot = X*0;
        for iComp = 1:K
          Ftmp = gm{it,K}.ComponentProportion(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
          F{iComp} = reshape(Ftmp,size(X));  
          Ftot = Ftot + F{iComp};
          Fplot = squeeze(sum(F{iComp},[1 2]))*ntot(it)*dvx*dvy*1e3;  
          plot(hca,zvec,Fplot)  
        end
        Fplot = squeeze(sum(Ftot,[1 2]))*ntot(it)*dvx*dvy*1e3;     
        plot(hca,zvec,Fplot)
        hca.XLabel.String = 'v_z (km/s)';
        hca.YLabel.String = 'f(v_z) (s/m^4)';
        irf_legend(hca,["Observed" "iComp = "+(1:K) "All comp."]',[0.98 0.98],'fontsize',fontsize)
        irf_legend(hca,["Observed","T = "+round(T_scalar{it,K})+" eV", "..."]',[0.02 0.98],'fontsize',fontsize)
      end
      hold(hca,'off')
    end
  
    if 1 % Plot component fits for one K
      K = 4;
      for it = 1:nt
        hca = h2(isub); isub = isub + 1;
        time = times(it);
        pdist = PD.tlim(time+0.5*0.150*[-1 1]);
        vdf = pdist.reduce('1D',[0 0 1],'nMC',nMC);
        plot(hca,vdf.depend{1},vdf.data,'linewidth',1.5)
        hold(hca,'on')
        mu = gm{it,K}.mu;
        Sigma = gm{it,K}.Sigma;
        gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm{it,K},[x0 y0 z0]),x,y,z);
        
        Ftot = X*0;
        for iComp = 1:K
          Ftmp = gm{it,K}.ComponentProportion(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
          F{iComp} = reshape(Ftmp,size(X));  
          Ftot = Ftot + F{iComp};
          Fplot = squeeze(sum(F{iComp},[1 2]))*ntot(it)*dvx*dvy*1e3;  
          plot(hca,zvec,Fplot)  
        end
        Fplot = squeeze(sum(Ftot,[1 2]))*ntot(it)*dvx*dvy*1e3;     
        plot(hca,zvec,Fplot)
        hca.XLabel.String = 'v_z (km/s)';
        hca.YLabel.String = 'f(v_z) (s/m^4)';
        irf_legend(hca,["Observed" "iComp = "+(1:K) "All comp."]',[0.98 0.98],'fontsize',fontsize)
        irf_legend(hca,["Observed","T = "+round(T_scalar{it,K})+" eV", "..."]',[0.02 0.98],'fontsize',fontsize)
      end
      hold(hca,'off')
    end
  
    if 0 % BIC
      hca = h2(isub); isub = isub + 1;
      bic = cellfun(@(x)x.BIC,gm);
      plot(hca,1:nGroupsMax,bic,'*')
      hca.XLim = [0.5 nGroupsMax+0.5];
      hca.XLabel.String = 'K';
      hca.YLabel.String = 'BIC';
      irf_legend(hca,"VDF"+(1:nt),[0.98 0.98],'fontsize',fontsize)
    end
    if 0 % AIC and BIC
      hca = h2(isub); isub = isub + 1;
      aic = cellfun(@(x)x.AIC,gm);
      bic = cellfun(@(x)x.BIC,gm);
      plot(hca,1:nGroupsMax,aic,'*')
      hca.XLim = [0.5 nGroupsMax+0.5];
      hca.XLabel.String = 'K';
      hca.YLabel.String = 'AIC';
      hold(hca,'on')
      plot(hca,1:nGroupsMax,bic,'s')
      hold(hca,'off')  
    end
  
    c_eval('axis(h2(?),''square'');',1:numel(h2))
    compact_panels(h2,0.05,0.05,1)
    drawnow
    if doPrint 
      cn.print(sprintf('gmm_vdf_df_id%04.0f_t0_%s_K234',iDF,time.utc('yyyymmdd_HHMMSS')))
    end
  end
  catch
    disp(sprintf('Skipping idf=%g due to error.',iDF))
  end
end

%% Some post-processing
db = table_ci(table_ci.id_df>0,:); % remove all empty values
ids = db.id_df;
for ids = 1:numel(ids)
  T_scalar_before =  db.T_scalar_before(ids);
  T_scalar_after  =  db.T_scalar_after(ids);
  T_max_before = cellfun(@(x)max(x),T_scalar_before);
  T_max_after  = cellfun(@(x)max(x),T_scalar_after);
  T_max_before(ids(id)) = T_max_before;
  T_max_after(ids(id)) = T_max_after;
end

%% Plot some general statistics :)
db = table_ci(table_ci.id_df>0,:); % remove all empty values
%  check if difference in temperature between components is a good
%  indicator

h = setup_subplots(3,2,'horizontal');
isub = 1;

hca = h(isub); isub = isub + 1;
histogram(hca,min(cat(1,db.T_min_after{:})'),0:50:2000)
hca.XLabel.String = 'min(T)';
hca.YLabel.String = 'Counts';

hca = h(isub); isub = isub + 1;
histogram(hca,max(cat(1,db.T_min_after{:})'),0:50:2000)
hca.XLabel.String = 'min(T)';
hca.YLabel.String = 'Counts';


hca = h(isub); isub = isub + 1;
T_min_plot = cat(1,db.T_min_after{:});
plot(hca,db.id_df,T_min_plot','-')
hca.XLabel.String = 'id';
hca.YLabel.String = 'min(T)';

if 0
hca = h(isub); isub = isub + 1;
histogram(hca,min(cat(1,db.T_min_after{:})'),0:50:2000)
hca.XLabel.String = 'min(T)';
hca.YLabel.String = 'Counts';
end

hca = h(isub); isub = isub + 1;
minT = min(cat(1,db.T_min_after{:})');
maxT = max(cat(1,db.T_min_after{:})');
toplot = minT./maxT;
histogram(hca,toplot,0:0.05:1)
hca.XLabel.String = 'min(T)/max(T)';
hca.YLabel.String = 'Counts';

hca = h(isub); isub = isub + 1;
minT = min(cat(1,db.T_min_after{:})');
maxT = max(cat(1,db.T_min_after{:})');
idx = kmeans([minT', maxT'],4);
idx = idx*0;

scatter(hca,minT,maxT,20,idx,'filled')
hca.XLabel.String = 'min(T)';
hca.YLabel.String = 'max(T)';
axis(hca,'square')
axis(hca,'equal')
hca.XLim(1) = 0;
hca.YLim(1) = 0;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';


hca = h(isub); isub = isub + 1;
maxTafter = min(cat(1,db.T_max_after{:})');
maxTbefore = max(cat(1,db.T_max_before{:})');
idx = kmeans([minT', maxT'],4);
idx = idx*0;

scatter(hca,maxTafter,maxTbefore,20,'k','filled')
hca.XLabel.String = 'max(T_{after})';
hca.YLabel.String = 'max(T_{before})';
axis(hca,'square')
axis(hca,'equal')
hca.XLim(1) = 0;
hca.YLim(1) = 0;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';

%%
hca = h(isub); isub = isub + 1;

T_min_plot = cat(1,db.T_min_before{:});
plot(hca,db.id_df,T_min_plot','-')
hca.XLabel.String = 'id';
hca.YLabel.String = 'min(T)';