% get tints
sep.get_tints; % tint_waves

% load additional data if needed
%c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
%c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',1:4);
%c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',1:4)
[Rpar,Rperp,thetaBR] = mms.mms4_displacement('gseB?','gseR?','plot',0);

% filter electric field
fmin = 10;
fmax = 0;
filter_order = 5;
c_eval('E?par = gseE?par.filt(fmin,fmax,[],filter_order).tlim(tint_waves);',1:4)
sampling_time_E = (E1par.time(2) - E1par.time(1)); % s

% set boundary parameters
max_time_delay = 3*1e-3; % s
min_abs_E = 5; % mV/m
max_time_p2p = 5*1e-3; % s
step_tint = 0.5*max_time_p2p; % s
length_tint = 0.5*max_time_p2p;
max_lags = ceil(max_time_delay/sampling_time_E);

% set up for loop for sub time intervals
n_tints = ceil((tint_waves.stop - tint_waves.start)/step_tint);
correlation_pairs = nchoosek(1:4,2);
n_correlation_pairs = size(correlation_pairs,1);
matrix_correlation = zeros(n_tints,n_correlation_pairs);
matrix_time_lags = zeros(n_tints,n_correlation_pairs);
matrix_time_ind = zeros(n_tints,n_correlation_pairs);
matrix_sc_distance_par = zeros(n_tints,n_correlation_pairs);
matrix_sc_distance_perp = zeros(n_tints,n_correlation_pairs);
tint_step_center = '';

%[ind1 ind2] = ind12(n,n_tints); % perhaps quicker to get the time indices like this first, instead of using tlim for each time interval below

fprintf('it = %4.0f/%4.0f\n',0,n_tints) % display progress
for i_tint = 1:n_tints
  if mod(i_tint,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],i_tint,n_tints); end % display progress
  if i_tint == n_tints
    tint_tmp = [tint_waves.start + (i_tint-1)*step_tint tint_waves.stop];
  else
    tint_tmp = tint_waves.start + (i_tint-1)*step_tint + [0 length_tint];
  end
  tint_tmp_center = tint_tmp.start + 0.5*(tint_tmp.stop-tint_tmp.start);
  tint_step_center(i_tint,:) = tint_tmp_center.utc;
  
  for i_pair = 1:n_correlation_pairs
    c_eval('E1 = E?par.tlim(tint_tmp).data;',correlation_pairs(i_pair,1));
    c_eval('E2 = E?par.tlim(tint_tmp).data;',correlation_pairs(i_pair,2));
    if numel(E1)>numel(E2) % make sure they have the same number of elements
      E1 = E1(1:numel(E2));
    elseif  numel(E2)>numel(E1)
      E2 = E2(1:numel(E1));
    end  
    if any([max(E1)<min_abs_E max(E1)<min_abs_E])
      continue;
    end
    [C,LAGS] = xcorr(E1,E2,max_lags);
    [max_C,ind_max] = max(C);
    matrix_time_ind(i_tint,correlation_pairs(i_pair,1),correlation_pairs(i_pair,2)) = ind_max;
    if ~isinteger(ind_max)
      ind_max = uint8(ind_max);
    end
    lag_at_max_C = LAGS(ind_max);
    matrix_correlation(i_tint,i_pair) = max_C;
    matrix_time_lags(i_tint,i_pair) = lag_at_max_C;
  end  
end

%% Make TSeries, for plotting and inspecting
tint_step_center = EpochTT(tint_step_center);
tsC = irf.ts_scalar(tint_step_center,matrix_correlation);
tsLag = irf.ts_scalar(tint_step_center,matrix_time_lags);
tsRpar = irf.ts_scalar(tint_step_center,Rpar.resample(tint_step_center).data);
tsRperp = irf.ts_scalar(tint_step_center,Rperp.resample(tint_step_center).data);
dt = matrix_time_lags*sampling_time_E;
tsVph = irf.ts_scalar(tint_step_center,tsRpar.data./dt);

lim_good_correlation = 4000; 
[ind_good_mean_correlation,good_mean_correlation] = find(mean(matrix_correlation,2)>lim_good_correlation);
[ind_good_correlation,good_correlation] = find(matrix_correlation>lim_good_correlation);
dt_good_correlation = dt(ind_good_mean_correlation,:);
mean_dt_good_correlation = mean(dt_good_correlation,1);
mean_dt_good_correlation = [0 mean_dt_good_correlation(1:3)];
dt_ref_mms1 = -mean_dt_good_correlation(1:4);

tsVph.data(ind_good_correlation) = NaN;

%% Plot
npanels = 8;
h = irf_plot(npanels);
str_correlation_pairs = cell(n_correlation_pairs,1);
for i_pair = 1:n_correlation_pairs; str_correlation_pairs{i_pair,1} = sprintf('%g%g',correlation_pairs(i_pair,1),correlation_pairs(i_pair,2)); end;
if 1 % E single sc
  hca = irf_panel('E 1sc');
  irf_plot(hca,E1par)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 1 % E 4 sc
  hca = irf_panel('E 4sc');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{E1par,E2par,E3par,E4par},'comp')
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 1 % E 4 sc time shifted
  hca = irf_panel('E 4sc dt');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{E1par,E2par,E3par,E4par},'comp','dt',dt_ref_mms1)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 1 % Correlation
  hca = irf_panel('Correlation');
  irf_plot(hca,tsC)
  hca.YLabel.String = {'C',''};
  irf_legend(hca,str_correlation_pairs,[1.01 0.999])
end
if 1 % Lag
  hca = irf_panel('Lags');
  irf_plot(hca,tsLag)
  hca.YLabel.String = {'Lag',''};
end
if 1 % dRpar
  hca = irf_panel('dR par');
  irf_plot(hca,tsRpar,'.')
  hca.YLabel.String = {'\Delta_{||} R','(km)'};
end
if 1 % dRperp
  hca = irf_panel('dR perp');
  irf_plot(hca,tsRperp,'.')
  hca.YLabel.String = {'\Delta_{\perp} R','(km)'};
end
if 1 % vph
  hca = irf_panel('vph');
  irf_plot(hca,tsVph,'*')
  hca.YLabel.String = {'v_{ph}','(km/s)'};
end

irf_zoom(h,'x',tint_waves)

c_eval('h(?).YLabel.Interpreter = ''tex'';',1:npanels)

