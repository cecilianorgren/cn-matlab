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
  %%
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

mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');

db_table_df = db_table_ff(db_table_ff.is_df==1,:);
nDF = numel(db_table_df.time);


iDFs = [9 128];
iDFs = 217:nDF;
%iDFs = 87;
iDFs = 17;23:nDF;
iDFs = 101:nDF;
iDFs = [1 6 9 10 13 21 25 35 38];
iDFs = 37;
iDFs = 57;
doPrint = 0;
doPlot = 1;
for iDF = iDFs%87%iDFs(1)
  %try
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
  c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
  iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;

  %% Prep data
  %PD_use = iPDist_counts; % iPDist
  PD_use = iPDist; % iPDist

  nMovMean = 7; % Number of distributions for the moving average. Use the same for finding the energy limt and the gmm

  %tsElow = PD_use.find_noise_energy_limit(5).movmean(15);
  tsElow = iPDist_counts.find_noise_energy_limit_counts(5,nMovMean);
  emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit
  PD = PD_use.mask({emask_mat});
  PD = PD.movmean(nMovMean);
  PD = PD(1:nMovMean:PD.length);
%  delete(hca)
  irf_plot(PD.omni.deflux.specrec); hca = gca; hca.YScale = 'log';
  
  %% Do the Gaussian Mixture Model
  tint_df = tDF + [-5 100]; % apply to two times, before and after DF
  times = PD.tlim(tint_df).time.start:nMovMean*0.150:PD.tlim(tint_df).time.stop;
  nMP = 100000; % Number of macroparticles  
  
  for Ks = [2 3]
    nGroupsMax = Ks;  % number of classes/groups for kmeans and gmm
    K = Ks;
  
    % Initial guess
    initial_guess = [0 -1000 -1000;...
                     0 -1000 +1000;...
                     0 +1000 -1000;...
                     0 +1000 +1000;...
                     0  0000  0000;...
                     0  0000  0000]*0.1;
  
    
 
    nt = times.length; 
    
    %T_mat = cell(nt,nGroupsMax);
    %T_scalar = cell(nt,nGroupsMax);
    clear gm
    tic
    for it = 1:nt
      if mod(it,10)==0; disp([it+"/" + nt]); end
      time = times(it);
      pdist = PD.tlim(time+0.5*0.150*[-1 1]);
      %nMP = nansum(round(pdist.data(~isnan(pdist.data)))); % total number of counts
      
      scpot = mean(scPot.tlim(time + 0.5*0.15*[-1 1]).data,1);
      scpot = irf.ts_scalar(time,scpot);
      B = mean(gseB.tlim(time + 0.5*0.15*[-1 1]).data,1); B = B/norm(B);
      
      for nGroups = nGroupsMax %1:nGroupsMax
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
            
        if it == 1
          S = 'randSample';
        else
          clear S
          S.mu = gm{it-1,K}.mu;
          S.Sigma = gm{it-1,K}.Sigma;
          S.ComponentProportion = gm{it-1,K}.ComponentProportion;        
          %S.mu = initial_guess(1:nGroups,:);
          %S.Sigma = repmat([500 10 10; 10 500 10; 10 10 500],[1 1 nGroups]);  
          %S.ComponentProportion = repmat(1,[1,nGroups]);
          if any(S.ComponentProportion==0); S = 'randSample'; end
        end
        
  
        %options = statset('Display','final','MaxIter',1500,'TolFun',1e-5);
        options = statset('MaxIter',150,'TolFun',1e-5);
  
  
        X = [MP.vx, MP.vy, MP.vz];
        %X(sqrt(sum(X.^2,2))>1000,:) = [];
        %gm{it,nGroups} = fitgmdist(X,nGroups,'Start',S,'SharedCovariance',false,'Options',options);
        gm{it,nGroups} = fitgmdist(X,nGroups,'Start',S,'SharedCovariance',false,'Options',options,'RegularizationValue',0.1);
  
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
    toc
    % Diagnose cold-ion quantities
    % Find coldest temperature for each it, K
    %T_min = cellfun(@(x)min(x),T_scalar);
    %T_max = cellfun(@(x)max(x),T_scalar);
  
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
  
    if doPlot
      %%
      fontsize = 13;
      vlim = 2500;
      nMC = 500;
    
      nvx = 101; nvy = 102; nvz = 103;
      xvec = linspace(-vlim,vlim,nvx); dvx = xvec(2)-xvec(1);
      yvec = linspace(-vlim,vlim,nvy); dvy = yvec(2)-yvec(1);
      zvec = linspace(-vlim,vlim,nvz); dvz = zvec(2)-zvec(1);
      [X,Y,Z] = ndgrid(xvec,yvec,zvec);
      XYZ = [X(:) Y(:) Z(:)];
  
      Ftot = zeros(size(X));
      %K = 4;
      clear fx fy fz vx vy vz Txx Tyy Tzz T_tens 
      T_tens = zeros(times.length,3,3,K);
      vx = zeros(times.length,K);
      vy = zeros(times.length,K);
      vz = zeros(times.length,K);
      n = zeros(times.length,K);
      for ii = 1:size(gm,1)
        [~,isort] = sort(gm{ii,K}.mu(:,3)); % vz
        [~,isort] = sort(squeeze(gm{ii,K}.Sigma(1,1,:)+gm{ii,K}.Sigma(2,2,:)+gm{ii,K}.Sigma(3,3,:))); % Ts
        mu = gm{ii,K}.mu(isort,:);
        Sigma = gm{ii,K}.Sigma(:,:,isort);
        compProp = gm{ii,K}.ComponentProportion(isort);
        for iComp = 1:K          
          Ftmp = compProp(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
          Ftmp = Ftmp;
          F{iComp} = reshape(Ftmp,size(X));  
          Ftot = Ftot + F{iComp};
          fx(ii,:,iComp) = sum(F{iComp},[2 3])*ntot(ii)*dvy*dvz*1e3;
          fy(ii,:,iComp) = sum(F{iComp},[1 3])*ntot(ii)*dvx*dvz*1e3;
          fz(ii,:,iComp) = sum(F{iComp},[1 2])*ntot(ii)*dvx*dvy*1e3;
  
          vx(ii,iComp) = mu(iComp,1);
          vy(ii,iComp) = mu(iComp,2);
          vz(ii,iComp) = mu(iComp,3);
  
          n(ii,iComp) = compProp(iComp)*ntot(ii) ;
  
          %vt2 = gm{it,nGroups}.Sigma*1e6; % (km/s)^2 -> (m/s)^2
          %cov_T_mat = units.mp*vt2/2/units.eV;
          %T_mat{it,nGroups} = cov_T_mat;
          T_tens(ii,:,:,:) = Sigma*1e6*units.mp/2/units.eV;
          %Tyy(ii,iComp) = gm{it,nGroups}.Sigma(2,2)*1e6*units.mp/2/units.eV;
          %Tzz(ii,iComp) = gm{it,nGroups}.Sigma(3,3)*1e6*units.mp/2/units.eV;
        end
      end
  %
  
      fx_tot =  sum(fx,3);
      fy_tot =  sum(fy,3);
      fz_tot =  sum(fz,3);
      
      clear tsFx tsFy tsFz tsT tsV tsN
      for iComp = 1:K
        tsFx{iComp} = PDist(times,fx(:,:,iComp),'1Dcart',xvec); tsFx{iComp}.units = 's/m^4';
        tsFy{iComp} = PDist(times,fy(:,:,iComp),'1Dcart',yvec); tsFy{iComp}.units = 's/m^4';
        tsFz{iComp} = PDist(times,fz(:,:,iComp),'1Dcart',zvec); tsFz{iComp}.units = 's/m^4';
        tsT{iComp} = irf.ts_tensor_xyz(times, T_tens(:,:,:,iComp));
        tsV{iComp} = irf.ts_vec_xyz(times, [vx(:,iComp),vy(:,iComp),vz(:,iComp)]);
        tsN{iComp} = irf.ts_scalar(times, n(:,iComp));
      end
      tsFx_tot =th PDist(times,fx_tot,'1Dcart',xvec);
      tsFy_tot = PDist(times,fy_tot,'1Dcart',yvec);
      tsFz_tot = PDist(times,fz_tot,'1Dcart',zvec);
      
      n_vec = cellfun(@(x) x.data, tsN, 'UniformOutput', false);
      tsN_tot = irf.ts_scalar(times,sum(cat(2,n_vec{:}),2))
  
  
      h = irf_plot(7);
  
      hca = irf_panel('B');
      irf_plot(hca,gseB)
      hca.YLabel.String = 'B (nT)';
  
      hca = irf_panel('fvx mms');
      %vdfx = PD.reduce('1D',[1 0 0]);
      %vdfx = PD.reduce('1D',[1 0 0]);
      vdf = PD.reduce('1D',[0 0 1]);
      irf_spectrogram(hca,vdf.specrec,'log')   
      hca.YLabel.String = 'v_z (km/s)';  
      irf_legend(hca,'MMS',[0.98 0.98],'k')
     
      hca = irf_panel('Fz tot');
      specrec = tsFz_tot.specrec('velocity');
      irf_spectrogram(hca,specrec,'log')    
      hca.YLabel.String = 'v_z (km/s)';
      irf_legend(hca,'GMM',[0.98 0.98],'k')
  
      if 1 % Tscalar of individual components      
        hca = irf_panel('T comp');
        ts = cellfun(@(x){x.xx+x.yy+x.zz},tsT);
        irf_plot(hca,ts,'comp')      
        hca.YLabel.String = 'T (eV)';       
        hca.YScale ='log';
        irf_zoom(hca,'y')
        irf_legend(hca,arrayfun(@(x) "Comp " + x,1:K,'UniformOutput',false)',[1.01 0.98])
      end
      if 1 % v of coldest (first) component
        hca = irf_panel('v 1');
        iComp = 1;
        tsVcold = irf.ts_vec_xyz(tsV{iComp}.time,[tsV{iComp}.x.data tsV{iComp}.y.data tsV{iComp}.z.data]);      
        irf_plot(hca,{tsVcold.x, tsVcold.y, tsVcold.z},'comp')
        hca.YLabel.String = sprintf('v_{%g} (km/s)',iComp);
        irf_legend(hca,{'v_x','v_y','v_z'},[0.98 0.98])
      end
      if 1 % v of second coldest (second) component
        hca = irf_panel('v 2');
        iComp = 2;
        tsVcold = irf.ts_vec_xyz(tsV{iComp}.time,[tsV{iComp}.x.data tsV{iComp}.y.data tsV{iComp}.z.data]);      
        irf_plot(hca,{tsVcold.x, tsVcold.y, tsVcold.z},'comp')
        hca.YLabel.String = sprintf('v_{%g} (km/s)',iComp);
        irf_legend(hca,{'v_x','v_y','v_z'},[0.98 0.98])
      end
      if 1 % n of individual components      
        hca = irf_panel('n comp');
        irf_plot(hca,tsN,'comp');      
        hold(hca,'on')
        irf_plot(hca,tsN_tot,'comp','k');
        hold(hca,'off')
        hca.YLabel.String = 'n (cm^{-3})';
        hca.YLabel.Interpreter = 'tex';
      end
      if 0 % vx of individual components      
        hca = irf_panel('vx comp');
        tsVx = cellfun(@(x){x.x},tsV);
        irf_plot(hca,tsVx,'comp')      
        hca.YLabel.String = 'v_x (km/s)';        
      end
      if 0 % vy of individual components      
        hca = irf_panel('vy comp');
        tsVy = cellfun(@(x){x.y},tsV);
        irf_plot(hca,tsVy,'comp')      
        hca.YLabel.String = 'v_y (km/s)';        
      end
      if 0 % vz of individual components      
        hca = irf_panel('vz comp');
        tsVz = cellfun(@(x){x.z},tsV);
        irf_plot(hca,tsVz,'comp')      
        hca.YLabel.String = 'v_z (km/s)';        
      end
  
  
      if 0 % f(vz) of individual components
        for iComp = 1:K
          hca = irf_panel(sprintf('Fz %g',iComp));
          specrec = tsFz{iComp}.specrec('velocity');
          irf_spectrogram(hca,specrec,'log')      
        end
      end
      
      hlinks = linkprop(h(2:3),{'CLim','YLim'});
  
      h(1).Title.String = sprintf('K = %g',K);
      colormap(flipdim(irf_colormap(hca,"spectral"),1))
      h(end).XTickLabelRotation = 0;
      irf_zoom(h,'x',[tsN{1}.time.start tsN{1}.time.stop]+[-5 5])
      irf_plot_axis_align
      c_eval('h(?).FontSize = 12;',1:numel(h))    
      hl = findobj(gcf,'type','line')
      c_eval('hl(?).LineWidth = 1;',1:numel(hl))
      c_eval('h(?).LineWidth = 1;',1:numel(h))
  
      drawnow
      cn.print(sprintf('gmm_iDF=%04.f_K=%g',iDF,K))
    end

  end
  %catch
  %  disp(sprintf('Skipping idf=%g due to error.',iDF)) 
  %end
end

%% Some post-processing
db = table_ci(table_ci.id_df>0,:); % remove all empty values
ids = db.id_df;
T_max_before = {};
T_max_after = {};
T_min_comp_before = {};
T_min_comp_after = {};
for id = 1:numel(ids)
  T_scalar_before =  db.T_scalar_before(id);
  T_scalar_after  =  db.T_scalar_after(id);
  T_max_before{ids(id),1} = cellfun(@(x)max(x),T_scalar_before{:});
  T_max_after{ids(id),1} = cellfun(@(x)max(x),T_scalar_after{:});
end

%% Some post-processing, find the coldest component, not coldest scalar
db = table_ci(table_ci.id_df>0,:); % remove all empty values
ids = db.id_df;
T_min_comp_before = cell(size(table_ci,1),1);
T_min_comp_after = cell(size(table_ci,1),1);
for id = 1:numel(ids)
  T_mat_before =  db.T_before(id);
  T_mat_after  =  db.T_after(id);
  maxK = size(T_mat_before{1},2);
  for K = 1:maxK
    diag_el_before = [];
    diag_el_after = [];
    T_mat_before_tmp = T_mat_before{1}{K};
    T_mat_after_tmp = T_mat_before{1}{K};
    for iK = 1:K
      diag_el_before = [diag_el_before diag(T_mat_before_tmp(:,:,iK))];
      diag_el_after = [diag_el_before diag(T_mat_after_tmp(:,:,iK))];
    end

    T_min_comp_before_tmp(K) = min(diag_el_before(:));
    T_min_comp_after_tmp(K) = min(diag_el_after(:));
    1;
  %T_max_before{ids(id),1} = cellfun(@(x)max(x),T_scalar_before{:});
  %T_max_after{ids(id),1} = cellfun(@(x)max(x),T_scalar_after{:});
  end
  T_min_comp_before(ids(id),1) = {T_min_comp_before_tmp};
  T_min_comp_after(ids(id),1) = {T_min_comp_after_tmp};
end

table_ci = addvars(table_ci,T_min_comp_before,T_min_comp_after);
%table_ci = removevars(table_ci,["T_min_comp_before","T_min_comp_after"]);
%% Some post-processing, remove K=5 which was there by mistake
%db = table_ci(table_ci.id_df>0,:); % remove all empty values
if 0
%for id = 1:numel(ids)
  for ivar = 5:size(table_ci,2)
    
    var = table_ci.(ivar);
    %var = var{1};
    if iscell(var)
      table_ci.(ivar) = cellfun(@(x) x(1:end-1),var,'UniformOutput',false);
    elseif isnumeric(var)
      table_ci(it,ivar) = var(1:4);
    end
  end
%end
end
%% Plot some general statistics :)
db = table_ci(table_ci.id_df>0,:); % remove all empty values
%  check if difference in temperature between components is a good
%  indicator

h = setup_subplots(3,3,'horizontal');
isub = 1;

hca = h(isub); isub = isub + 1;
histogram(hca,min(cat(1,db.T_min_after{:})'),0:100:2000)
hca.XLabel.String = 'min_K(min_i((T_{after}))';
hca.YLabel.String = 'Counts';
yyaxis(hca,'right')
hca = gca;

hca = h(isub); isub = isub + 1;
histogram(hca,max(cat(1,db.T_min_after{:})'),0:50:2000)
hca.XLabel.String = 'max_K(min_i(T_{after}))';
hca.YLabel.String = 'Counts';

hca = h(isub); isub = isub + 1;
T_min_plot = cat(1,db.T_max_before{:});
plot(hca,db.id_df,T_min_plot','-')
hca.XLabel.String = 'id';
hca.YLabel.String = 'max_i(T_{before})';

hca = h(isub); isub = isub + 1;
T_min_plot = cat(1,db.T_min_after{:});
plot(hca,db.id_df,T_min_plot','-')
hca.XLabel.String = 'id';
hca.YLabel.String = 'min(T_{after})';


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
minT = min(cat(1,db.T_min_before{:})');
maxT = max(cat(1,db.T_min_before{:})');
idx = kmeans([minT', maxT'],4);
idx = idx*0+1;
scatter(hca,minT,maxT,20,idx,'filled')
hold(hca,'on')
maxlim = max([hca.XLim, hca.YLim]);
plot(hca,[0 maxlim],[0 maxlim],'k')
hold(hca,'off')
hca.XLabel.String = 'min(T_{before})';
hca.YLabel.String = 'max(T_{before})';
axis(hca,'equal')
axis(hca,'square')
hca.XLim(1) = 0;
hca.YLim(1) = 0;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';
colormap(hca,[0 0 0])

hca = h(isub); isub = isub + 1;
minT = min(cat(1,db.T_min_after{:})');
minT_comp = min(cat(1,db.T_min_comp_after{:})');
scatter(hca,minT,minT_comp,20,idx,'filled')
hold(hca,'on')
maxlim = max([hca.XLim, hca.YLim]);
plot(hca,[0 maxlim],[0 maxlim],'k')
hold(hca,'off')
hca.XLabel.String = 'min(T_{after})';
hca.YLabel.String = 'min(T_{after}^{comp})';
axis(hca,'equal')
axis(hca,'square')
hca.XLim(1) = 0;
hca.YLim(1) = 0;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';
colormap(hca,[0 0 0])

hca = h(isub); isub = isub + 1;
minT = min(cat(1,db.T_min_after{:})');
maxT = max(cat(1,db.T_max_after{:})');
minT_comp = min(cat(1,db.T_min_comp_after{:})');
%maxT = max(cat(1,db.T_max_after{:})');
idx = kmeans([minT', maxT'],4);
idx = idx*0+1;
scatter(hca,minT,maxT,20,idx,'filled')
hold(hca,'on')
scatter(hca,minT_comp,maxT,20,idx+2,'filled')
maxlim = max([hca.XLim, hca.YLim]);
plot(hca,[0 maxlim],[0 maxlim],'k')
hold(hca,'off')
hca.XLabel.String = 'min(T_{after})';
hca.YLabel.String = 'max(T_{after})';
axis(hca,'equal')
axis(hca,'square')
hca.XLim(1) = 0;
hca.YLim(1) = 0;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';
colormap(hca,[0 0 0])


hca = h(isub); isub = isub + 1;
maxTafter = max(cat(1,db.T_max_after{:})');
maxTbefore = max(cat(1,db.T_max_before{:})');
idx = kmeans([minT', maxT'],4);
idx = idx*0;
scatter(hca,maxTafter,maxTbefore,20,'k','filled')
hold(hca,'on')
maxlim = max([hca.XLim, hca.YLim]);
plot(hca,[0 maxlim],[0 maxlim],'k')
hold(hca,'off')
hca.XLabel.String = 'max(T_{after})';
hca.YLabel.String = 'max(T_{before})';
axis(hca,'equal')
axis(hca,'square')
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