%% Spacecraft, time interval, location
savePath = '/Users/cecilia/Research/dipolarization_fronts/initial_search';

irf.log('critical')
units = irf_units;
ic_order = [2 3 4 1];
tint_all = irf.tint('2018-01-01T00:00:00.00Z/2019-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% criteria for location
Rcrit = {[-Inf -8],[-15 15],[-10 10]};

% save('/Volumes/Fountain/Data/PIC/df_save_all_2017.mat','df_save_all','Rcrit','criteria')
%% Define criteria for dipolarization fronts
% Bz = B0 + B1*tanh(t/dt)
dt = 1; % s, initial guess
%dBmin = 5; % nT

% moving window
w_step = 0.1; % s, window step
w_size = 8; % s, window size

% selection criteria for fit parameters and correlation
% Two main options for parameter definitions
% B0 - middle value of front
% B1 - half change of B across front
% dt - two half times of front
% Bfit = @(t,B0,B1,dt) B0 + B1*tanh(t/dt);
criteria = {[0 4],[-3 15],[5 Inf],[0.9 Inf]}; % 2*dt, B0-B1, B1*2, xcorr
% c - B before front, requirement > 0
% a - change in B across the front
% dt - duration of front, two half lengths/widths/times
% setting requirement for B value ahead of front is easier for latter
% version, in the first version, it depends on both B0 and B1
%Bfit = @(t,B0,dB,dt) (c+a/2) + (a/2)*tanh(t/(dt/2));
%criteria = {[0 2],[0 20],[3 Inf],[0.9 Inf]}; % dt, B0, B1, xcorr



%% Run through files
% initialize empty table
df_save_all = [];
nfiles = numel(files);

for ifile = 1:nfiles
  fileId = strsplit(files(ifile).name,'_'); fileId = fileId{5};
  tint = [files(ifile).start files(ifile).stop];
  strdisp = sprintf('#%g/%g, Time: %s - %s.',ifile,nfiles,tint(1).utc,tint(2).utc);
  disp(strdisp)    
  
  %% Load spacecraft location and check location criteria
  R = mms.get_data('R_gsm',tint);
  
  if isfield(R,'gsmR1') && size(R.gsmR1,2) == 4
    c_eval('gsmR? = irf.ts_vec_xyz(R.time,R.gsmR?(:,2:4));',1:4); % dfg_srvy_l2pre
  elseif isfield(R,'gsmR1')
    c_eval('gsmR? = irf.ts_vec_xyz(R.time,R.gsmR?);',1:4); % mec
  else
    disp('  No position data. Skipping time interval.')
    continue;
  end
  if     not(isempty(gsmR1)), gsmR = gsmR1;
  elseif not(isempty(gsmR2)), gsmR = gsmR2;
  elseif not(isempty(gsmR3)), gsmR = gsmR3;
  elseif not(isempty(gsmR4)), gsmR = gsmR4;  
  else    
    disp('  No position data. Skipping time interval.')
    continue
  end
  Rx = gsmR(1).x.data/(units.RE*1e-3);
  Ry = gsmR(1).y.data/(units.RE*1e-3);
  Rz = gsmR(1).z.data/(units.RE*1e-3);
  if Rx > Rcrit{1}(1) && ...
     Rx < Rcrit{1}(2) && ...
     Ry > Rcrit{2}(1) && ...
     Ry < Rcrit{2}(2) && ...
     Rz > Rcrit{3}(1) && ...
     Rz < Rcrit{3}(2)   
    strdisp = sprintf('  R = [%.2f, %.2f, %.2f] RE.',Rx,Ry,Rz);
    disp(strdisp)
  else    
    strdisp = sprintf('  R = [%.2f, %.2f, %.2f] RE. Aborting.',Rx,Ry,Rz);
    disp(strdisp)    
    continue;
  end  
    
  %% Load data: gsmB*  
  ok = 0;
  iic = 0;
  while not(ok) && iic < 5
    iic = iic + 1; 
    ic = ic_order(iic);
    c_eval('gsmB = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic); 
    if isempty(gsmB)
      ok = 0;
    else
      ok = 1;
    end
  end
  if ok == 0
    disp('  Found no B-field data. Aborting.')    
    continue
  end
  disp(sprintf(' . Using MMS %g',ic))
  
  %% Look for DFs
  % run through time interval
  df_data = df_detection(gsmB,dt,w_step,w_size);
  % apply criteria
  if not(isempty(df_data))
    df_save = df_selection(df_data,criteria);
  else
    df_save = [];
  end
  
  strdisp = sprintf('  Found %g potential DF(s).',numel(df_save));
  disp(strdisp)    
  
  %% If a df was found, load some more data and save in table
  if not(isempty(df_save))
    tint_load = tint + [1 -1];
    c_eval('gseE = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint_load);',ic);
    gsmE = c_coord_trans('GSE','GSM',gseE);
    c_eval('ne = mms.get_data(''Ne_fpi_brst_l2'',tint_load,?);',ic);
    c_eval('ni = mms.get_data(''Ni_fpi_brst_l2'',tint_load,?);',ic);
    c_eval('gseVe = mms.get_data(''Ve_gse_fpi_brst_l2'',tint_load,?);',ic)
    c_eval('gseVi = mms.get_data(''Vi_gse_fpi_brst_l2'',tint_load,?);',ic);    
    gsmVe = c_coord_trans('GSE','GSM',gseVe);
    gsmVi = c_coord_trans('GSE','GSM',gseVi);
    
    for isave = 1:numel(df_save)
      df_save(isave).ic = ic;
      tint_save = [df_save(isave).t1 df_save(isave).t2];      
      df_save(isave).data.fileId = fileId;      
      df_save(isave).data.gsmR = gsmR.tlim(tint_save);
      df_save(isave).data.gsmB = gsmB.tlim(tint_save);
      if not(isempty(gsmE)), df_save(isave).data.gsmE = gsmE.tlim(tint_save); else, df_save(isave).data.gsmE = []; end
      if not(isempty(gsmVe)), df_save(isave).data.gsmVe = gsmVe.tlim(tint_save); else, df_save(isave).data.gsmVe = []; end
      if not(isempty(gsmVi)), df_save(isave).data.gsmVi = gsmVi.tlim(tint_save); else, df_save(isave).data.gsmVi = []; end
      if not(isempty(ni)), df_save(isave).data.ni = ni.tlim(tint_save); else, df_save(isave).data.ni = []; end
      if not(isempty(ne)), df_save(isave).data.ne = ne.tlim(tint_save); else, df_save(isave).data.ne = []; end
            
      df_save(isave).R = [Rx Ry Rz];
    end
    df_save_all = cat(2,df_save_all,df_save);
  end  
  %% Plot results
  if 1 % plot results
    %%    
    npanels = 4;
    h = irf_plot(npanels);
    
    if 1 % B
      hca = irf_panel('B');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmB.x,gsmB.y,gsmB.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'B_{GSM} (nT)'};
    end
    if not(isempty(gsmE)) % E
      hca = irf_panel('E');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmE.x,gsmE.y,gsmE.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'E_{GSM} (mV/m)'};
    end
    if not(isempty(ne)) % ne
      hca = irf_panel('n');
      set(hca,'ColorOrder',mms_colors('12'))
      irf_plot(hca,{ne},'comp')
      set(hca,'ColorOrder',mms_colors('12'))
      irf_legend(hca,{'n_e'},[0.98 0.98])
      hca.YLabel.String = {'n (cm^{-3})'};
    end
    if not(isempty(gsmVi)) % Vi
      hca = irf_panel('Vi');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmVi.x,gsmVi.y,gsmVi.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'V_{i,GSM} (km/s)'};
    end
    
    c_eval('irf_pl_mark(h(?),[df_save(!).t1 df_save(!).t2]);',1:npanels,1:numel(df_save))   
    %c_eval('irf_pl_mark(h(?),df_save(!).t1 + 0.5*(df_save(!).t2-df_save(!).t1),''color'',[0 0 0]);',1:npanels,1:numel(df_save))    
    irf_zoom(h,'x',tint)
    
    if not(isempty(df_save))
      str_print = sprintf('df_save_%s_nDF=%02.0f',df_save(1).data.fileId,numel(df_save));      
    else
      str_print = sprintf('df_save_%s_nDF=%02.0f',fileId,0);
    end
    cn.print(str_print,'path','/Users/cecilia/Research/dipolarization_fronts/initial_search')
  end
  
end
disp('Done for now.')

%% Plot results 
nrows = 3;
ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;

df = df_save_all; % for ease of notation
ndf = numel(df);

all_opt_best = [df.opt_best];
opt_best_B0 = all_opt_best(1,:); 
opt_best_B1 = all_opt_best(2,:); 
opt_best_dt = all_opt_best(3,:); 
opt_xcorr = [df.opt_xcorr];
all_R = cat(1,df.R);

if 1 % plot all Bx
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0)
  hold(hca,'on')
  for idf = 1:ndf        
    tt = df(idf).data.gsmB.time-df(idf).data.gsmB.time(1);
    bb = df(idf).data.gsmB.x.data;
    plot(hca,tt,bb)
  end    
  hold(hca,'off')
  hca.XTick = 0:6;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B_{x,GSM}';
  %legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1 % plot all By
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0)
  hold(hca,'on')
  for idf = 1:ndf        
    tt = df(idf).data.gsmB.time-df(idf).data.gsmB.time(1);
    bb = df(idf).data.gsmB.y.data;
    plot(hca,tt,bb)
  end    
  hold(hca,'off')
  hca.XTick = 0:6;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B_{y,GSM}';
  %legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1 % plot all Bx
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0)
  hold(hca,'on')
  for idf = 1:ndf        
    tt = df(idf).data.gsmB.time-df(idf).data.gsmB.time(1);
    bb = df(idf).data.gsmB.z.data;
    plot(hca,tt,bb)
  end    
  hold(hca,'off')
  hca.XTick = 0:6;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B_{z,GSM}';
  %legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1 % plot all n
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0)
  hold(hca,'on')
  for idf = 1:ndf        
    if not(isempty(df(idf).data.ne))
      tt = df(idf).data.ne.time-df(idf).data.ne.time(1);
      bb = df(idf).data.ne.data;
      plot(hca,tt,bb)
    else
      continue
    end
  end    
  hold(hca,'off')
  hca.YLabel.String = 'n_{e}';
  %legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1 % Fit parameters, B1 vs B0
  hca = h(isub); isub = isub + 1;
  plot(hca,opt_best_B0,opt_best_B1,'*')
  hold(hca,'on')
  plot(hca,[0 max([hca.XLim hca.YLim])],[0 max([hca.XLim hca.YLim])],'k')
  hold(hca,'off')
  hca.YLabel.String = 'B_1';
  hca.XLabel.String = 'B_0';
  hca.Title.String = 'B = B0 + B1*tanh(t/dt)';
end
if 1 % Fit parameters B1 vs dt
  hca = h(isub); isub = isub + 1;
  plot(hca,opt_best_dt,opt_best_B1,'*')
  hca.YLabel.String = 'B_1';
  hca.XLabel.String = 'dt';
  hca.Title.String = 'B = B0 + B1*tanh(t/dt)';
end
if 0 % Position, X, sqrt(Y2+Z2)
  hca = h(isub); isub = isub + 1;
  plot(hca,df(idf).R(1),sqrt(df(idf).R(2).^2 + df(idf).R(3).^2),'*')
  hca.XLabel.String = 'R_x';
  hca.YLabel.String = '(R_y^2 + R_z^2)^{1/2}';  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % Position, X, Y
  hca = h(isub); isub = isub + 1;
  plot(hca,all_R(:,1),all_R(:,2),'*')
  hca.XLabel.String = 'R_x';
  hca.YLabel.String = 'R_y';  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % Position, X, Z
  hca = h(isub); isub = isub + 1;
  plot(hca,all_R(:,1),all_R(:,3),'*')
  hca.XLabel.String = 'R_x';
  hca.YLabel.String = 'R_z';  
  hca.YLabel.Interpreter = 'tex';
end



