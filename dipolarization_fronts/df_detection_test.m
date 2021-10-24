%% Check all files I have on in database locally
tint_all = irf.tint('2015-01-01T00:00:00.00Z/2022-01-01T00:00:00.00Z');
allfiles = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

%%
ic = 1;
tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('gsmE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gsm_brst_l2'',tint);',ic);
c_eval('ne? = mms.get_data(''Ne0_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)

%% Test df_detection()
dt = 1; % s
w_step = 0.1; % s, window step
w_size = 6; % s
dBmin = 5; % nT
tic
df_data = df_detection(gsmB1,dBmin,dt,w_step,w_size);
toc
criteria = {[0 3],[2 Inf],[0.9 Inf]}; % dt, B, xcorr
df_save = df_selection(df_data,criteria);

%%
nrows = 6;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

t = [df_data.it];
all_opt_best = [df_data.opt_best];
opt_best_B0 = all_opt_best(1,:); 
opt_best_B1 = all_opt_best(2,:); 
opt_best_dt = all_opt_best(3,:); 
opt_xcorr = [df_data.opt_xcorr];

% make some kind of selection
selection = zeros(size(t));
ind = opt_best_dt<3 & opt_best_B1>2 & opt_xcorr>0.9;
selection(ind) = 1;

% find consecutive sequences and find larges correlation within each of
% them
[B, N, BI] = RunLength(selection);
% find all with slection==1
idf = [];
all_seq_sel = find(B==1);
for isec = 1:numel(B)
  if B(isec) == 1
    part_data = opt_xcorr(BI(isec):(BI(isec)+N(isec)-1));
    [maxval,i_max_corr] = max(part_data);
    idf(end+1) = BI(isec)+i_max_corr-1;
  else continue, end
end

% Save all data in some readable way
df_save = df_data(idf);

if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,t,[df_data.reg_RMS],...
           t,[df_data.opt_RMS])
  hca.YLabel.String = 'C_{rms}';
  legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,t,[df_data.reg_xcorr],...
           t,[df_data.opt_xcorr])
  hca.YLabel.String = 'C_{xcorr}';
  legend(hca,{'lin reg','opt'},'location','best','box','off')
end
if 1
  hca = h(isub); isub = isub + 1;  
  plot(hca,t,opt_best_B0)
  hca.YLim = [0 20];
  hca.YLabel.String = 'B_{0}^{opt}';
end
if 1
  hca = h(isub); isub = isub + 1;  
  plot(hca,t,opt_best_B1)
  hca.YLim = [-10 10];
  hca.YLabel.String = 'B_{1}^{opt}';
end
if 1
  hca = h(isub); isub = isub + 1;  
  plot(hca,t,opt_best_dt)
  hca.YLim = [0 5];
  hca.YLabel.String = 'dt^{opt}';
end
if 1 % selection
  hca = h(isub); isub = isub + 1;  
  plot(hca,t,selection,'-',idf,ones(size(idf)),'*')
  hca.YLim = [0 1.1];
  hca.YLabel.String = 'selection';
end
if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,t,[df_data.theta])
  hleg = legend(hca,{'B0','B1'});
  hleg.Title.String = sprintf('B, = B0 + B1*tanh(t/%.2f)',dt);
  hleg.Box = 'off';
end

irf_plot_axis_align
hlinks = linkprop(h,{'XLim'});
hlinks.Targets(1).XLim = [0 t(end)];

