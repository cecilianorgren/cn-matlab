%% Load data
ic = 1:4;
tint = irf.tint('2017-07-06T08:16:03.00Z/2017-07-06T08:18:13.00Z');

% Load datastore
mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   

c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
disp('Loading spacecraft position...')
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)


%% Correlation
% 0.01 s duration perhaps
dt_sampling_original = gseE1.time(2)-gseE1.time(1);
dt_window = 0.01;
window_lengths_t = dt_window*[1:2:10];
n_step_levels = numel(window_lengths_t);
window_lengths_t_steps = fix(window_lengths_t/dt_sampling_original);
step_window = window_lengths_t_steps;

%tintZoom  = irf.tint('2017-07-06T08:16:35.00Z',10);
tintZoom = irf.tint('2017-07-06T08:16:37.00Z',3);
c_eval('E = gseE?par.tlim(tintZoom);',ic)
times_E = E.time;
n_times_E = E.time.length;

for ilevel = 1:n_step_levels
  % number of windows
  n_windows = n_times_E/window_lengths_t_steps(ilevel)
  for iwindow = 1:n_windows
    
    timeline = gseE1.tlim(tint_vicinity); timeline = tint_vicinity(1):0.5*dt_sampling_original:tint_vicinity(2);
    c_eval('totE? = gseE?.tlim(tint_vicinity).resample(timeline);')
    c_eval('totE? = totE?.resample(timeline);')
    c_eval('E? = gseE?par.tlim(tint_vicinity).resample(timeline);')
    c_eval('E? = E?.resample(timeline);')
    c_eval('R? = gseR?.resample(timeline).tlim(tint);')

    dt_sampling = E1.time(2)-E1.time(1);
    dt = zeros(4,1);
    C = ones(4,1);
    for ic = 2:4
      c_eval('[tmpC,lags] = xcorr(E1.data,E?.data,''coeff'');',ic)  
      i_shift = find(abs(tmpC) == max(abs(tmpC)));
      C(ic) = tmpC(i_shift);
      di = -lags(i_shift);
      dt(ic) = di*dt_sampling;
    end

    c_eval('matR? = [R?.time.epochUnix R?.data];',1:4)
    v_xcorr = irf_4_v(matR1,matR2,matR3,matR4,dt + E1(1).time.epochUnix); %v_xcorr = v_xcorr(2:4);
    v_direction = irf_norm(v_xcorr);
    v_amplitude = sqrt(sum(v_xcorr.^2));

  end
end


