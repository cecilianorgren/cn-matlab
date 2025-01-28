% wave_analysis_timing_redo

%% Load time intervals
% fprintf(fid,[data_format '\n'],...
%   tint_vicinity(1).utc,tint_vicinity(2).utc,... % 1 2
%   tint_esw(1).utc,tint_esw(2).utc,...           % 3 4
%   t_center.utc,...                              % 5
%   v_xcorr,...                                   % 6 7 8
%   v_amplitude,...                               % 9
%   esw.phi_max,...                               % 10 11 12 13
%   esw.peaktopeak,...                            % 14
%   esw.C...                                      % 15 16 17 18
% );
data_format = '%s %s %s %s %s %.1f %.1f %.1f %.1f %4.0f %4.0f %4.0f %4.0f %2.1f %.3f %.3f %.3f %.3f'; 
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 

fid = fopen([matlabPath 'esw_properties.txt'],'r');
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);
end

doPlot = 0;
doWrite = 1;
doOverwrite = 1;
%doAppend = 0;
fileWrite = 'esw_properties_redo_dt.txt';
dtmult = 1;
% dt: c_eval('gseEdt? = irf_integrate(E?,tint_esw(1)+1*dt(?));');  % dtmult = 1 or 0
if doOverwrite % wipe file
  fid = fopen([matlabPath fileWrite],'w');  
end
%%
% Data is appended to [matlabPath 'esw_properties_redo.txt'], so remember 
% to remove data from that file if this is redone, or if the script stops
% in the middle for some reason. If not data will be duplicated.
redoesw = 1:numel(esw_data{end});
for iesw = 1:numel(redoesw)  
  ii = redoesw(iesw);
  
  %% Get time
  tint_vicinity = EpochTT([esw_data{1}{ii}; esw_data{2}{ii}]);
  tint_esw = EpochTT([esw_data{3}{ii}; esw_data{4}{ii}]);
  t_center = EpochTT(esw_data{5}{ii});
  
  %% Get ESW v, and potential
  dt_sampling_original = gseE1.time(2)-gseE1.time(1);
  timeline = gseE1.tlim(tint_vicinity); 
  timeline = tint_vicinity(1):0.5*dt_sampling_original:tint_vicinity(2);
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

  c_eval('tsV?_timing = irf.ts_vec_xyz(totE?.time,repmat(v_xcorr,totE?.length,1));')

  % integrate only Epar and multiply with v_amplitude
  direction_sign = gseB1.resample(t_center).norm.dot(v_direction).data;
  c_eval('gseEdt? = irf_integrate(E?,tint_esw(1)+dtmult*dt(?));');  
  c_eval('gsePhi? = gseEdt?*v_amplitude*direction_sign;')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')

  c_eval('v_trapping? = sqrt(2*units.e*gsePhi?.abs/units.me)*1e-3;')

  l = fit_gaussian(E1,v_amplitude);

  esw = struct;
  esw.tint_vicinity = tint_vicinity;
  esw.tint_esw = tint_esw;
  esw.t_center = t_center;
  esw.phi_max = [max(findpeaks(gsePhi1.tlim(tint_esw).data)) max(findpeaks(gsePhi2.tlim(tint_esw).data)) max(findpeaks(gsePhi3.tlim(tint_esw).data)) max(findpeaks(gsePhi4.tlim(tint_esw).data))];
  esw.v_trapping_max = sqrt(2*units.e*esw.phi_max/units.me)*1e-3;
  esw.peaktopeak = l*2;
  esw.debye_length = [Ld1.resample(t_center).data Ld2.resample(t_center).data Ld3.resample(t_center).data Ld4.resample(t_center).data];
  esw.C = C;
  %[mean(esw.phi_max) std(esw.phi_max)]
  %esw.phi_max
  
  if 0 % plot
    h = irf_plot(2);
    colors = pic_colors('matlab');
    
    hca = irf_panel('E par');
    irf_plot(hca,'E?','comp','dt',0*dt)
    hca.YLabel.String = 'E par (mV/m)';
    
    hca = irf_panel('phi');
    set(hca,'ColorOrder',colors)
    irf_plot(hca,'gsePhi?','comp','dt',0*dt)
    for ii=1:4; hpl = irf_pl_mark(hca,tint_esw(1)+1*dt(ii)); hpl.Color = colors(ii,:); end
    hca.YLabel.String = '\phi (V)';
    pause
  end
    
  
  %% Write data to file
  if doWrite
    fid = fopen([matlabPath fileWrite],'a+');
    fprintf(fid,[data_format '\n'],...
      tint_vicinity(1).utc,tint_vicinity(2).utc,... % 1 2
      tint_esw(1).utc,tint_esw(2).utc,...           % 3 4
      t_center.utc,...                              % 5
      v_xcorr,...                                   % 6 7 8
      v_amplitude,...                               % 9
      esw.phi_max,...                               % 10 11 12 13
      esw.peaktopeak,...                            % 14
      esw.C...                                      % 15 16 17 18
    );
  end
end
