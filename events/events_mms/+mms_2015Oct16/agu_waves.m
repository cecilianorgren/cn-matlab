%% Waves

tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:05.00Z'); % magnetosphere-magnetosheath-magnetosphere

%c_eval('plE? = irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(L).data dslE?brst.dot(M).data dslE?brst.dot(N).data]);')
c_eval('plE? = dslE?brst;')

h = irf_plot(4);

ic  = 4;
% Magnetic field
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.7]);

% Electric field
hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[1.01 0.7]);

% E
hca=irf_panel('E wavelet'); %irf_plot(hca,dslE1brst);
if 1 % wavelet  
  c_eval('irf_spectrogram(hca,wavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('E pfft');
  c_eval('irf_spectrogram(hca,pfftE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end
% B
  hca=irf_panel('B wavelet'); %irf_plot(hca,dmpaB1scm);
if 1 % wavelet
  c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('B pfft');
  c_eval('irf_spectrogram(hca,pfftB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end

irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h(1:6),'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end
%% Waves

tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:05.00Z'); % magnetosphere-magnetosheath-magnetosphere

%c_eval('plE? = irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(L).data dslE?brst.dot(M).data dslE?brst.dot(N).data]);')
c_eval('plE? = dslE?brst;')

h = irf_plot(8);

ic  = 1;
% Magnetic field
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.7]);

% Electric field
hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[1.01 0.7]);

% E
hca=irf_panel('E wavelet'); irf_plot(hca,dslE1brst);
if 0 % wavelet  
  c_eval('irf_spectrogram(hca,wavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('E pfft');
  c_eval('irf_spectrogram(hca,pfftE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end
% B
  hca=irf_panel('B wavelet'); irf_plot(hca,dmpaB1scm);
if 0 % wavelet
  c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
elseif 0 % pfft
  hca=irf_panel('B pfft');
  c_eval('irf_spectrogram(hca,pfftB?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?brst,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 -1]; 
  hca.YTick = [1e2 1e3 1e4];
  hca.YLim=[10 4080];
end

irf_zoom(h(1:4),'x',tint)
%irf_pl_mark(h(1:2),tint)
hca = irf_panel('delete for space'); %delete(hca)    
grid(hca,'off')    
irf_timeaxis(h(4))
    
if 1 %plot lh
  hca = irf_panel('lh phi')
  tint_lh = irf.tint('2015-10-16T10:33:24.00Z',1.5); 
  E = dslE3brst.tlim(tint_lh);
  acB = dmpaB3scm.tlim(tint_lh);
  dcB = dmpaB3.tlim(tint_lh);
  c_eval('n_loc = mean(ne_psd?.tlim(tint_lh).data);',sc)
  c_eval('Te_loc = mean(Te?perp.tlim(tint_lh).data);',sc)

  angles = 0:2:359;
  f_highpass = 0.5;

  % Propagation direction
  [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir({dcB,acB},E,angles,f_highpass);
  i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
  direction=x(i_dir,:);

  % Velocity and density
  mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3.mean

  v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
  phiB_scale = B0*1e-18/mu0/e/n; % km/s
  v = v_approx*linspace(0.1,2.5,30);

  [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
  i_v=find(corr_v==min(corr_v));
  velocity=v(i_v);
  tsBz = irf.ts_scalar(irf_time(Bz(:,1),'epoch>epochunix'),Bz(:,2));
  tsIntEdt = irf.ts_scalar(irf_time(intEdt(:,1),'epoch>epochunix'),intEdt(:,1+i_dir));
  irf_plot(hca,{tsBz*phiB_scale,tsIntEdt*velocity},'comp')
  %irf_legend(hca,{'B_z*scale','\int Edt*v'},[0.98 0.95])
  irf_legend(hca,{'\phi_{B}','\phi_{E}'},[0.98 0.95])
  hca.YLabel.String = '\Phi [V]'; 
  irf_legend(hca,{['\omega > ' num2str(flim,'%.2f')]},[0.02 0.1])
  irf_zoom(hca,'x',tint_lh)
  %irf_zoom(hca,'y',[-300 300])    
  axtop = add_length_on_top(hca,velocity,0.5);
  axtop.XLabel.String = 'Length [km]';    
  % add B labels on right
  yratio = phi_B(100,2)/Bz(100,2);            
  axtop.YLim = hca.YLim/yratio;
  axtop.YLabel.String = '\delta B_{||} [nT]'; 
  axtop.YTickLabelMode = 'auto';
  axtop.YTickMode = 'auto';
  % mark different physical lengths
  irf_plot_axis_align
  irf_plot_zoomin_lines_between_panels(h(4),h(6))
  
    
  % Compare Ek and Bz amplitudes
  maxEk = max(dEk(:,2));
  maxEn = max(dEn(:,2));
  maxBz = max(Bz(:,2));
  disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
end   

hca = irf_panel('delete for space 2'); %delete(hca)    
grid(hca,'off')    
irf_timeaxis(h(6))

tint_waves = irf.tint('2015-10-16T10:33:30.30Z',0.1); 
c_eval('Epar = irf.ts_scalar(dslE?brst.time,dslE?brst.dot(dmpaB1brst.resample(dslE?brst)).data)',ic)
velocity = 800;
phiEpar = irf_integrate(Epar.tlim(tint_waves))*velocity;
hca = irf_panel('Epar phi');
irf_plot(hca,phiEpar);

irf_plot_axis_align
irf_plot_zoomin_lines_between_panels(h(4),h(8))
hca = irf_panel('delete for space'); delete(hca) 
hca = irf_panel('delete for space 2'); delete(hca) 
irf_zoom(h,'x',irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z'))
irf_zoom(h(1:6),'y')
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end

%%
fig = figure('Position',[0.1300    0.1100    500    150],'units','normalized');
hca = irf_panel('lh phi');
tint_lh = irf.tint('2015-10-16T10:33:24.00Z',1.5); 
  E = dslE3brst.tlim(tint_lh);
  acB = dmpaB3scm.tlim(tint_lh);
  dcB = dmpaB3.tlim(tint_lh);
  c_eval('n_loc = mean(ne_psd?.tlim(tint_lh).data);',sc)
  c_eval('Te_loc = mean(Te?perp.tlim(tint_lh).data);',sc)

  angles = 0:2:359;
  f_highpass = 0.5;

  % Propagation direction
  [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir({dcB,acB},E,angles,f_highpass);
  i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
  direction=x(i_dir,:);

  % Velocity and density
  mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3.mean

  v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
  phiB_scale = B0*1e-18/mu0/e/n; % km/s
  v = v_approx*linspace(0.1,2.5,30);

  [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
  i_v=find(corr_v==min(corr_v));
  velocity=v(i_v);
  tsBz = irf.ts_scalar(irf_time(Bz(:,1),'epoch>epochunix'),Bz(:,2));
  tsIntEdt = irf.ts_scalar(irf_time(intEdt(:,1),'epoch>epochunix'),intEdt(:,1+i_dir));
  set(hca,'ColorOrder',mms_colors('xy'))
  irf_plot(hca,{tsBz*phiB_scale,tsIntEdt*velocity},'comp')
  %irf_legend(hca,{'B_z*scale','\int Edt*v'},[0.98 0.95])
  set(hca,'ColorOrder',mms_colors('xy'))
  irf_legend(hca,{'\phi_{B}','\phi_{E}'},[0.02 0.95])
  irf_legend(hca,{'c'},[0.98 0.95],'color','black','fontsize',14)
  hca.YLabel.String = '\phi [V]'; 
  
  %irf_legend(hca,{['\omega > ' num2str(flim,'%.2f')]},[0.02 0.1])
  irf_legend(hca,{['f > ' num2str(flim,'%.2f')]},[0.02 0.1])
  irf_zoom(hca,'x',tint_lh)
  %irf_zoom(hca,'y',[-300 300])    
  axtop = add_length_on_top(hca,velocity,0.5);
  axtop.XLabel.String = 'Length [km]';    
  % add B labels on right
  yratio = phi_B(100,2)/Bz(100,2);            
  axtop.YLim = hca.YLim/yratio;
  axtop.YLabel.String = '\delta B_{||} [nT]'; 
  axtop.YTickLabelMode = 'auto';
  axtop.YTickMode = 'auto';
  axtop.Position(4)=axtop.Position(4)*0.75;
  hca.Position(4)=hca.Position(4)*0.75;
  axtop.Position(2)=axtop.Position(2)*1.15;
  hca.Position(2)=hca.Position(2)*1.15;
  axtop.Position(3)=axtop.Position(3)*0.9;
  hca.Position(3)=hca.Position(3)*0.9;
  axtop.FontSize = 14;
  hca.FontSize = 14;
  