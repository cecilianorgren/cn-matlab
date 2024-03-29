
times = get_time(4,'epochtt'); times = times.sort;
tint_vicinity = times([1 4]);
tint_esw = times([2 3]);
t_center = tint_esw(1) + 0.5*(tint_esw(2)-tint_esw(1));
irf_pl_mark(gca,tint_esw);

tint_str = [irf_time(tint_esw(1),'epochtt>utc_yyyymmdd_HHMMSSmmm') '-' irf_time(tint_esw(2),'epochtt>utc_HHMMSSmmm')];

%% Get ESW v, and potential
dt_sampling_original = gseE1.time(2)-gseE1.time(1);
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

c_eval('tsV?_timing = irf.ts_vec_xyz(totE?.time,repmat(v_xcorr,totE?.length,1));')

if 0 % integrate total E and dit with v
  c_eval('gseEdt? = irf_integrate(totE?,tint_esw(1));');
  c_eval('gseEdt? = irf.ts_vec_xyz(gseEdt?.time,gseEdt?.data);')
  phi_filt = 3;
  c_eval('gsePhi? = gseEdt?.dot(tsV?_timing); gsePhi?_filt = gsePhi?.filt(phi_filt,0,[],10);')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
else % integrate only Epar and multiply with v_amplitude
  c_eval('gseEdt? = irf_integrate(E?,tint_esw(1));');  
  c_eval('gsePhi? = gseEdt?*v_amplitude*double(sign(v_direction(1)));')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
end
  
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


%% Plot 
figure(10)
npanels = 7;
h = irf_plot(npanels);

hca = irf_panel('E');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{E1,E2,E3,E4},'comp')
hca.YLabel.String = 'E (mV/m)';
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  

hca = irf_panel('E dt');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{E1,E2,E3,E4},'comp','dt',dt)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('12341'))
irf_legend(hca,{['\Delta t = [ ' num2str(dt(1)*1e3)],num2str(dt(2)*1e3),num2str(dt(3)*1e3),num2str(dt(4)*1e3),'] ms'},[0.95 0.99],'fontsize',12);
set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('v = %.1f x [%.2f %.2f %.2f] km/s',v_amplitude,v_direction)},[0.95 0.05],'fontsize',12);
set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('C_{1>1234} = [%.2f %.2f %.2f %.2f]',C)},[0.01 0.99],'fontsize',12);


if 1 % Ex
  hca = irf_panel('Ex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.tlim(tint_vicinity).x,gseE2.tlim(tint_vicinity).x,gseE3.tlim(tint_vicinity).x,gseE4.tlim(tint_vicinity).x},'comp')
  hca.YLabel.String = {'E_x','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ey
  hca = irf_panel('Ey');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.tlim(tint_vicinity).y,gseE2.tlim(tint_vicinity).y,gseE3.tlim(tint_vicinity).y,gseE4.tlim(tint_vicinity).y},'comp')
  hca.YLabel.String = {'E_y','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ez
  hca = irf_panel('Ez');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.tlim(tint_vicinity).z,gseE2.tlim(tint_vicinity).z,gseE3.tlim(tint_vicinity).z,gseE4.tlim(tint_vicinity).z},'comp')
  hca.YLabel.String = {'E_z','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  
end

%hca = irf_panel('E dt by eye');
%irf_plot(hca,{E1,E2,E3,E4},'comp','dt',[0.0000    0.0004    0.0002    0.0002])

%hca = irf_panel('int(E) dt');
%irf_plot(hca,{gseEdt1,gseEdt2,gseEdt3,gseEdt4},'comp','dt',dt)

if 1
  hca = irf_panel('phi dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gsePhi1,gsePhi2,gsePhi3,gsePhi4},'comp','dt',dt)
  hca.YLabel.String = {'\Phi','(V)'};
end
if 0
  hca = irf_panel('phi dt detrend');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gsePhi1_detrend,gsePhi2_detrend,gsePhi3_detrend,gsePhi4_detrend},'comp','dt',dt)
  hca.YLabel.String = {'\Phi','(V)'};
end


hca = irf_panel('v trapping');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{v_trapping1*1e-3,v_trapping2*1e-3,v_trapping3*1e-3,v_trapping4*1e-3},'comp','dt',dt)
hca.YLabel.String = {'v_{trapping}','(10^3 km/s)'};

add_length_on_top(h(1),v_amplitude,1)
irf_zoom(h(:),'x',tint_vicinity)
irf_zoom(h(:),'y')
irf_plot_axis_align

cn.print(sprintf('esw_%s',tint_str),'path',eventPath)

c_eval('hmark(?,!) = irf_pl_mark(h(?),tint_esw(!),''k'');',1:npanels,1:2)


%% Write (append) data to file

fid = fopen([matlabPath 'esw_properties.txt'],'a+');
data_format = '%s %s %s %s %s %.1f %.1f %.1f %.1f %4.0f %4.0f %4.0f %4.0f %2.1f %.3f %.3f %.3f %.3f'; 

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
fid = fclose(fid);
print_str = ...
  sprintf([data_format],...
  tint_vicinity(1).utc,tint_vicinity(2).utc,... % 1 2
  tint_esw(1).utc,tint_esw(2).utc,...           % 3 4
  t_center.utc,...                              % 5
  v_xcorr,...                                   % 6 7 8
  v_amplitude,...                               % 9
  esw.phi_max,...                               % 10 11 12 13
  esw.peaktopeak...                             % 14
);
fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
