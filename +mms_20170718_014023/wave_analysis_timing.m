times = get_time(4,'epochtt'); times = times.sort;
tint_vicinity = times([1 4]);
tint_esw = times([2 3]);
t_center = tint_esw(1) + 0.5*(tint_esw(2)-tint_esw(1));

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
for ic = 2:4
  c_eval('[C,lags] = xcorr(E1.data,E?.data,''coeff'');',ic)
  i_shift = find(abs(C) == max(abs(C)));
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
irf_legend(hca,{sprintf('v = %.1f x [%.2f %.2f %.2f] km/s',v_amplitude,v_direction)},[0.95 0.01],'fontsize',12);

if 1 % Ex
  hca = irf_panel('Ex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.tlim(tint_vicinity).x,gseE2.tlim(tint_vicinity).x,gseE3.tlim(tint_vicinity).x,gseE4.tlim(tint_vicinity).x},'comp')
  hca.YLabel.String = {'E_x','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ex
  hca = irf_panel('Ey');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.tlim(tint_vicinity).y,gseE2.tlim(tint_vicinity).y,gseE3.tlim(tint_vicinity).y,gseE4.tlim(tint_vicinity).y},'comp')
  hca.YLabel.String = {'E_y','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ex
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


%% Write data to file

fid = fopen([matlabPath 'esw_properties.txt'],'a+');
data_format = '%s %s %s %s %s %.1f %.1f %.1f %.1f %4.0f %4.0f %4.0f %4.0f %2.1f'; 
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f'; 
fprintf(fid,[data_format '\n'],...
  tint_vicinity(1).utc,tint_vicinity(2).utc,...
  tint_esw(1).utc,tint_esw(2).utc,...
  t_center.utc,...
  v_xcorr,...
  v_amplitude,...
  esw.phi_max,...
  esw.peaktopeak...
)
fid = fclose(fid);


%% Load esw data
fid = fopen([matlabPath 'esw_properties.txt'],'r');
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);
end
  
fid = fopen([matlabPath 'esw_properties_redo.txt'],'r');
  
tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
tsVtrap = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{10}/units.me)*1e-3);
c_eval('tsVtrap? = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{9+?}/units.me)*1e-3);',1:4)

%%
tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
tsVtrap = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{10}/units.me)*1e-3);

npanels = 5;

h = irf_plot(npanels);
isub = 0;

if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{tsPhi});
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
end

if 0 % v phase
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{tsVphpar},'comp');
  hca.YLabel.String = {'v_{ph}','(km/s)'};  
end

if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase and trapping velocity');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{tsVphpar},'comp');
  hold(hca,'on')
  vmin = tsVphpar-tsVtrap;
  vmax = tsVphpar+tsVtrap;
  irf_patch(hca,{vmin,vmax})
  hca.YLim = sort([max(vmax.data) (min(vmin.data))]);
  hold(hca,'off')
  hca.YLabel.String = {'v_{ph}','(km/s)'};  
end

irf_zoom(h,'y')
irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));

%% Plot, including f proj and v phi and vtrap
ic = 1;
npanels = 10;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
tint_zoom = irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z');

vmin = tsVphpar-tsVtrap;
vmax = tsVphpar+tsVtrap;
  
  
if 0 % B abs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B abs');
  irf_plot(hca,gseB1.abs);
  hca.YLabel.String = 'B (nT)';
end
if 1 % B GSE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % iDEF omni
  isub = isub + 1;
  hca = irf_panel('iDEF');  
  [hout,hcb] = irf_spectrogram(hca,iDist.convertto('s^3/m^6').omni.specrec,'log');
%   hold(hca,'on')
%   lineScpot = irf_plot(hca,scPot1,'k');
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % iPDist pa 64
  isub = isub + 1;
  hca = irf_panel('i PA e64 deflux lowe');  
  eint = [100 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 1 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('iLine');
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVi},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_i (km/s)'; 
end
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.convertto('s^3/m^6').omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot_lim,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};   
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVe},'comp')
  set(hca,'ColorOrder',mms_colors('122'))
  irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_e (km/s)'; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{tsPhi});
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
end
if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  
  irf_spectrogram(hca,ef1D.specrec('velocity_1D'));
  hold(hca,'on')
  vmin = tsVphpar-tsVtrap;
  vmax = tsVphpar+tsVtrap;
  
  set(hca,'ColorOrder',mms_colors('122'))
  irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  hca.YLim = sort([max(vmax.data) (min(vmin.data))]);
  hold(hca,'off')
  hca.YLabel.String = {'v','(km/s)'};  
  set(hca,'ColorOrder',mms_colors('122'))
  irf_legend(hca,{'v_{ph}','v_{trap}'},[0.98 0.9],'fontsize',12);
  
  hsurf = findobj(hca.Children,'Type','Surface');
  hsurf.FaceAlpha = 0.2;
  hline = findobj(hca.Children,'Type','Line');
  c_eval('hline(?).LineWidth = 1.5;',1:numel(hline))
end

h(5).CLim = [-7 -3];
irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(4).CLim = [-35 -28]+12
colormap('jet');