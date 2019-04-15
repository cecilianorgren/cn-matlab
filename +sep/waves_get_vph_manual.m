
times = get_time(8,'epochtt'); %times = times.sort;
tint_vicinity = times([1 2]);
tint_esw = times([3 4]);
t_ref = times(5:8); % where to start put tref for integration ?
t_center = tint_esw(1) + 0.5*(tint_esw(2)-tint_esw(1));
hmark_tmp = irf_pl_mark(gca,tint_esw,'r'); hmark_tmp.FaceAlpha = 0.5;
hmark_tmp = irf_pl_mark(gca,t_ref.epochUnix,'k');

tint_str = [irf_time(tint_esw(1),'epochtt>utc_yyyymmdd_HHMMSSmmm') '-' irf_time(tint_esw(2),'epochtt>utc_HHMMSSmmm')];

sc_observing = zeros(1,4);
sc_observing_input = input('Which spacecraft observe the wave? [1234]:','s');
if isempty(sc_observing_input)
  sc_observing_input = '1234';
end
sc_observing = zeros(1,4);
for iisc = 1:4
  if strfind(sc_observing_input,num2str(iisc))
    sc_observing(iisc) = 1;
  end
end
       
%% Get ESW v, and potential
dt_sampling_original = gseE1.time(2)-gseE1.time(1);
timeline = gseE1.tlim(tint_vicinity); timeline = tint_vicinity(1):0.5*dt_sampling_original:tint_vicinity(2);
timeline_esw = gseE1.tlim(tint_esw); timeline_esw = tint_esw(1):0.5*dt_sampling_original:tint_esw(2);
c_eval('totE? = gseE?.tlim(tint_vicinity).resample(timeline);')
c_eval('totE? = totE?.resample(timeline);')
ffilt = 10;
c_eval('E? = gseE?par.tlim(tint_vicinity).resample(timeline).filt(ffilt,0,[],5);')
c_eval('E? = E?.resample(timeline);')
c_eval('R? = gseR?.resample(timeline_esw).tlim(tint_esw);')
c_eval('R?! = R?-R!;',1:4,1:4)
c_eval('B?! = 0.5*(gseB?.resample(timeline_esw)+gseB!.resample(timeline_esw));',1:4,1:4)
B1234 = 0.25*(gseB1.resample(timeline_esw)+gseB2.resample(timeline_esw)+gseB3.resample(timeline_esw)+gseB4.resample(timeline_esw));
avB1234 = mean(B1234.data);
c_eval('R?!par = mean(R!?.dot(B!?.norm).data); R_par_all(!,?) = R?!par;',1:4,1:4)

dt_sampling = E1.time(2)-E1.time(1);
dt_all = nan(4,4);
C_all = nan(4,4);
maxlag = round((tint_esw.stop-tint_esw.start)/dt_sampling/2);
iC = 1;
for ic1 = 1:4
  for ic2 = (ic1+1):4
    iC = iC + 1;
    c_eval('[tmpC,lags] = xcorr(E?.data,E!.data,''coeff'',maxlag);',ic1,ic2)  
    i_shift = find((tmpC) == max((tmpC)));    
    C_all(ic1,ic2) = tmpC(i_shift);
    %C_all(ic2,ic1) = C_all(ic1,ic2);
    di = -lags(i_shift);
    dt_all(ic1,ic2) = di*dt_sampling; 
    %dt_all(ic2,ic1) = dt_all(ic1,ic2);
    fprintf('ic1 = %g, ic2 = %g, C = %.2f, n_shift = %g, dt = %g \n',ic1,ic2,tmpC(i_shift),i_shift,di*dt_sampling)
  end
end
dt = dt_all(1,1:4); dt(1) = 0;% 11,12,13,14
v_all = -R_par_all./dt_all;
v_all(abs(v_all)==Inf) = NaN;
v_all_ = v_all; v_all_(:,sc_observing==0) = NaN; v_all_(sc_observing==0,:) = NaN;
v_all_(abs(v_all_)==Inf) = NaN;
v_ccorr = nanmean(v_all_(:));
Clim = 0.7;
v_all_high_C = v_all; v_all_high_C(C_all<Clim) = NaN;

c_eval('matR? = [R?.time.epochUnix R?.data];',1:4)
v_xcorr = irf_4_v(matR1,matR2,matR3,matR4,dt + E1(1).time.epochUnix); %v_xcorr = v_xcorr(2:4);
v_xcorr_par = dot(v_xcorr,avB1234/norm(avB1234));
v_direction = irf_norm(v_xcorr);
v_amplitude = sqrt(sum(v_xcorr.^2));

c_eval('tsV?_timing = irf.ts_vec_xyz(totE?.time,repmat(v_xcorr,totE?.length,1));')

if 0 % integrate total E and dit with v
  c_eval('gseEdt? = irf_integrate(totE?,tint_esw(1));');
  c_eval('gseEdt? = irf.ts_vec_xyz(gseEdt?.time,gseEdt?.data);')
  phi_filt = 3;
  c_eval('gsePhi? = gseEdt?.dot(tsV?_timing); gsePhi?_filt = gsePhi?.filt(phi_filt,0,[],10);')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
elseif 0 % integrate only Epar and multiply with v_xcorr_par
  c_eval('gseEdt? = irf_integrate(E?,t_ref(?));');  
  c_eval('gsePhi? = gseEdt?*abs(v_xcorr_par)*double(sign(v_direction(1)));')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
else % integrate only Epar and multiply with v_ccorr
  c_eval('gseEdt? = irf_integrate(E?,t_ref(?));');  
  c_eval('gsePhi? = gseEdt?*(1*v_ccorr);')
  c_eval('gsePhi?_detrend = gsePhi?; gsePhi?_detrend.data = detrend(gsePhi?_detrend.data,''linear'');')
end
  
c_eval('v_trapping? = sqrt(2*units.e*gsePhi?.abs/units.me)*1e-3;')

c_eval('l? = fit_gaussian(E?,v_amplitude);')

esw = struct;
esw.tint_vicinity = tint_vicinity;
esw.tint_esw = tint_esw;
esw.t_center = t_center;
%esw.debye_length = [Ld1.resample(t_center).data Ld2.resample(t_center).data Ld3.resample(t_center).data Ld4.resample(t_center).data];
esw.C = [C_all(1,2:4) C_all(2,3:4) C_all(3,4)];
esw.dt = [dt_all(1,2:4) dt_all(2,3:4) dt_all(3,4)];
esw.Rpar = [R_par_all(1,2:4) R_par_all(2,3:4) R_par_all(3,4)];
esw.v = [v_all(1,2:4) v_all(2,3:4) v_all(3,4)];
esw.v_xcorr = v_amplitude*v_xcorr;
esw.phi_max = [max(findpeaks(gsePhi1.tlim(tint_esw).data)) max(findpeaks(gsePhi2.tlim(tint_esw).data)) max(findpeaks(gsePhi3.tlim(tint_esw).data)) max(findpeaks(gsePhi4.tlim(tint_esw).data))];
esw.v_trapping_max = sqrt(2*units.e*esw.phi_max/units.me)*1e-3;
esw.peaktopeak = [l1 l2 l3 l4]*2;

%% Plot 
figure(10)
npanels = 7;
h = irf_plot(npanels);

hca = irf_panel('E');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{E1,E2,E3,E4},'comp')
hca.YLabel.String = 'E (mV/m)';
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.01 0.9],'fontsize',12);  

set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('f > %.0f Hz',ffilt)},[0.98 0.99],'fontsize',12);

set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('<B> = [%.1f %.1f %.1f] nT',avB1234)},[0.01 0.99],'fontsize',12);
set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('<v_{ccorr,par}> = %.0f km/s',v_ccorr)},[0.01 0.84],'fontsize',12);

set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('v_{xcorr} = [%.0f %.0f %.0f] km/s',v_xcorr)},[0.98 0.3],'fontsize',12);
set(hca,'ColorOrder',mms_colors('1'))
irf_legend(hca,{sprintf('v_{xcorr,par} = %.0f km/s',v_xcorr_par)},[0.98 0.1],'fontsize',12);



hca = irf_panel('E dt');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{E1,E2,E3,E4},'comp','dt',dt)
hca.YLabel.String = {'E','(mV/m)'};
%set(hca,'ColorOrder',mms_colors('12341'))
%irf_legend(hca,{['\Delta t = [ ' num2str(dt(1)*1e3)],num2str(dt(2)*1e3),num2str(dt(3)*1e3),num2str(dt(4)*1e3),'] ms'},[0.95 0.99],'fontsize',12);
%set(hca,'ColorOrder',mms_colors('1'))
%irf_legend(hca,{sprintf('v = %.1f x [%.2f %.2f %.2f] km/s',v_amplitude,v_direction)},[0.95 0.05],'fontsize',12);
%set(hca,'ColorOrder',mms_colors('1'))
%irf_legend(hca,{sprintf('C_{123>1234} = [%.2f %.2f %.2f %.2f %.2f %.2f %.2f]',C)},[0.01 0.99],'fontsize',12);


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

cn.print(sprintf('esw_%s',tint_str),'path',pathFigures)

c_eval('hmark(?,!) = irf_pl_mark(h(?),tint_esw(!),''k'');',1:npanels,1:2)


%% Write (append) data to file

fid = fopen([pathData sprintf('/esw_prop_%s.txt',tint_str)],'w');
%data_format = '%s %s %s %s %s %.1f %.1f %.1f %.1f %4.0f %4.0f %4.0f %4.0f %2.1f %.3f %.3f %.3f %.3f'; 

fprintf(fid,[data_format_write  '\n'],...
  tint_vicinity(1).utc,tint_vicinity(2).utc,... % 1 2
  tint_esw(1).utc,tint_esw(2).utc,...           % 3 4
  t_center.utc,...  
  t_ref(1).utc,t_ref(2).utc,t_ref(3).utc,t_ref(4).utc,...% 5
  avB1234,...
  sc_observing,...
  [C_all(1,2:4) C_all(2,3:4) C_all(3,4)],... 
  [dt_all(1,2:4) dt_all(2,3:4) dt_all(3,4)],...
  [R_par_all(1,2:4) R_par_all(2,3:4) R_par_all(3,4)],...
  [v_all(1,2:4) v_all(2,3:4) v_all(3,4)],...        %
  v_ccorr,...
  v_xcorr,...
  v_xcorr_par,...                                   % 6 7 8
  esw.phi_max,...                               % 10 11 12 13
  esw.peaktopeak...                            % 14
);
fid = fclose(fid);

% print_str = ...
%   sprintf([data_format],...
%   tint_vicinity(1).utc,tint_vicinity(2).utc,... % 1 2
%   tint_esw(1).utc,tint_esw(2).utc,...           % 3 4
%   t_center.utc,...                              % 5
%   v_xcorr,...                                   % 6 7 8
%   v_amplitude,...                               % 9
%   esw.phi_max,...                               % 10 11 12 13
%   esw.peaktopeak...                             % 14
% );
%fprintf('printed to file %sesw_properties.txt: %s /n',pathData,print_str)
