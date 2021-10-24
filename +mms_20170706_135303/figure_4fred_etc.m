units = irf_units;
% Load datastore
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   

% MMS id and time interval
mms_id = 1;
tint = irf.tint('2017-07-06T13:54:05.520Z/2017-07-06T13:54:05.640Z');
c_eval('eDist = ePDist?;',mms_id)
c_eval('eDist_bgremoved = eDist_nobg?.tlim(tint);',mms_id)

% Load EH properties
% observed/measured properties (konrad)
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_vtrap_all = sqrt(2*units.e*obs_potential/units.me);
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);

% EDI parameters
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

% Remove background electrons
nPhoto = 00;
nSecond = 20;
if 0
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint, ...
             'Nphotoe_art', nPhoto, 'nSecondary', nSecond, 'ZeroNaN', 0);
elseif 0
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint,...
             'ZeroNaN', 0);  
end
%% Make reduced distribution
strTint = [irf_time(tint(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vint = [-Inf Inf];
vg = (-100:2:100)*1e3;
c_eval('eDist = ePDist?.tlim(tint);',mms_id)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 40;
nMC = 500;
tic; ef1D = eDist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

tic; ef1D_bgremoved = eDist_bgremoved.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2_bgremoved = eDist_bgremoved.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

% Make pitch angle spectrograms
ePitch = eDist.pitchangles(dmpaB1.resample(eDist),12);
ePitch_bgremoved = eDist_bgremoved.pitchangles(dmpaB1.resample(eDist_bgremoved),12); 
%%
for it = 1:4
  %%
  figure(21)
  %it = 3;
  nrows = 2;
  ncols = 5;
  npanels = ncols*nrows;
  for ipanel = 1:npanels
    h(ipanel) = subplot(nrows,ncols,ipanel);
  end
  isub = 1;

  if 1 % par
    hca = h(isub); isub = isub + 1;
    plot(hca,ef1D(it).depend{1}*1e-3,ef1D(it).data);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f (s/m^4)';
  end
  if 1 % parperp1
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp1(it).plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  end
  if 1 % par1perp2
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp2(it).plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % perp1perp2
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_perp1perp2(it).plot_plane(hca);
    hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % pitch angle spectrogram
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ePitch(it).plot_pad_polar(hca,'scpot',scpot);
    hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end

  if 1 % par, bg removed
    hca = h(isub); isub = isub + 1;
    plot(hca,ef1D_bgremoved(it).depend{1}*1e-3,ef1D_bgremoved(it).data);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f (s/m^4)';
  end
  if 1 % parperp1, bg removed
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp1_bgremoved(it).plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';
  end
  if 1 % par1perp2, bg removed
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp2_bgremoved(it).plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % perp1perp2, bg removed
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_perp1perp2_bgremoved(it).plot_plane(hca);
    hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % pitch angle spectrogram, bg removed
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ePitch_bgremoved(it).plot_pad_polar(hca,'scpot',scpot);
    %hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
    %hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  
  h(1).Title.String = {ef2D_parperp1(it).time.utc,'no background removed'};
  h(6).Title.String = {'background removed'};
  for ipanel = [1 6]
    h(ipanel).CLim = [-13 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg([1 end])*1e-3*0.5;
    h(ipanel).YLim = [0 4e-3];
  end  
  links_2d = linkprop(h([2:4 7:9]),{'XLim','YLim','CLim'});
  for ipanel = [2:4 7:9]
    h(ipanel).CLim = [-15 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg([1 end])*1e-3*1;
    h(ipanel).YLim = vg([1 end])*1e-3*1;
  end
  colormap(pic_colors('candy4'))
  cn.print(sprintf('f1D_2D_bgrem_it=%g_nPhoto=%g_nSecond=%g_clim',it,nPhoto,nSecond))
end

% Run instability anaysis
%[f_ins,params_ins] = mms_20170706_135303.get_f0(2);
%[f0,params_0] = mms_20170706_135303.get_f0(1);

%% 4 fred, vph/trapping range, EDI range, average flux from model, and f0 of model

if 0
  %%
colors = mms_colors('matlab');
v_fpi = mean(ef1D_bgremoved.depend{1},1);
ind0 = find(v_fpi==0);
doAverage = 1;
if doAverage
  f_fpi = mean(ef1D.data,1);
  f_fpi_bgrem = nanmean(ef1D_bgremoved.data,1);
  fpi_utc = ef1D.time([1 end]).utc;
  fpi_utc = sprintf('%s - %s',fpi_utc(1,12:23),fpi_utc(2,12:23));
else
  f_fpi = ef1D.data;
  f_fpi_bgrem = ef1D_bgremoved.data;
  fpi_utc = ef1D.time.utc;
end

f_fpi(:,ind0) = [];
f_fpi_bgrem(:,ind0) = [];
v_fpi(:,ind0) = [];

figure(113)
nrows = 1;
ncols = 1;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 0 % f fpi
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi,'-');
  hca.YLabel.String = {'f (s^1/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
  
  if doAverage
    str_lines = fpi_utc;
  else
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
  end
  irf_legend(hca,str_lines,[0.99 0.99])
  set(hca,'ColorOrder',zeros(10,3))
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % f fpi, rem bg
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi_bgrem,'-');
  hca.YLabel.String = {'f (s/m^4)'};
  hca.XLabel.String = {'v (10^3 km/s)'};  
  
  if doAverage
    str_lines = sprintf('FPI time interval:\n%s',fpi_utc);
  else
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
  end
  h_fpi_str = irf_legend(hca,str_lines,[0.99 0.99]);
  set(hca,'ColorOrder',zeros(10,3))
  h_fpi_str.HorizontalAlignment = 'left';
  h_fpi_str.Position(1) = 0.62;
  h_fpi_str.Position(2) = 0.95;
  h_fpi_str.BackgroundColor = [1 1 1];
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % <f> model
  hold(hca,'on')
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_vec*v_scale*1e-3,mean(Fabel_obs,1),'-','color',colors(2,:));
  hca.YLabel.String = {'f_e (s/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
      %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  %hca.XLim = 1.7*[-30 30];
  %hca.YLim = [0 3.5]*1e-3;
  hold(hca,'off')
end
if 1 % f0
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [f0,params_0] = mms_20170706_135303.get_f0(1);
    plot(hca,v_plot*v_scale,f0(v_plot*1e3,params_0.n,params_0.vd,params_0.vt),'color',colors(3,:))
    if 1
      str_info = {'f_0 input parameters:'
                ['T_{0}= [' sprintf('%g  ',params_0.T) '] eV'];...
                ['n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cm^{-3}'];...
                ['v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];...
                };
    elseif 0
      str_info = ['f_0 input parameters:\n ',...
                'T_{0}= [' sprintf('%g  ',params_0.T) '] eV\n ',...
                'n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cc\n ',...
                'v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];
    else
      str_info = sprintf('fpi time interval:\n%s',fpi_utc);
    end
    h_f0 = irf_legend(hca,str_info,[0.90,0.35]);
    c_eval('h_f0(?).HorizontalAlignment = ''left'';',1:4)
    c_eval('h_f0(?).Position(1) = 0.62;',1:4)
    %c_eval('h_f0(?).Position(2) = 0.7;',1:4)
  
    c_eval('h_f0(?).BackgroundColor = [1 1 1];',1:4)
    c_eval('h_f0(?).BackgroundColor = ''none'';',1:4)
    hold(hca,'off')
end
if 1 % v_edi
  hold(hca,'on')
  if 1 % plot dashed lines
    h_edi = plot(hca,[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim*0.20,'k--',-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--');
    %plot(hca,-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--')
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vph patch
  hold(hca,'on')
  all_vph = [mean(obs_velocity)-std(obs_velocity);...
             mean(obs_velocity);... 
             mean(obs_velocity)+std(obs_velocity)]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,all_vph*1e-3,hca.YLim,'k-.')
  else % plot patch                    
    hpatch = patch(hca,1e-3*[all_vph(1,1),all_vph(end,1) all_vph(end,1) all_vph(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.7 0.9]);
    hpatch.FaceAlpha = 0.4;
    hpatch.EdgeAlpha = 0.4;
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vtrap patch
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,plot_vtrap_minus*1e-3,hca.YLim,'k-.',plot_vtrap_plus*1e-3,hca.YLim,'k-.')
  else % plot patch
    hpatch = patch(hca,1e-3*[plot_vtrap_plus(1,1),plot_vtrap_plus(end,1) plot_vtrap_plus(end,1) plot_vtrap_plus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.8 0.8]);
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '--';

    hpatch = patch(hca,1e-3*[plot_vtrap_minus(1,1),plot_vtrap_minus(end,1) plot_vtrap_minus(end,1) plot_vtrap_minus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.8 0.8]);
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '--';
  end
  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
yloc_errobar = 2.3*1e-3;
if 0 % vph errorbar for std
  hold(hca,'on')  
  %plot(hca,mean(obs_velocity)*1e-3*[1 1],hca.YLim)
  errorbar(hca,mean(obs_velocity)*1e-3,yloc_errobar,std(obs_velocity)*1e-3,'horizontal')
  hold(hca,'off')
end
if 0 % vtrap errorbar for std
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  
  %plot(hca,nanmean(vtrap_plus)*1e-3*[1 1],hca.YLim,'-')
  errorbar(hca,nanmean(vtrap_plus)*1e-3,yloc_errobar,std(vtrap_plus,'omitnan')*1e-3,'horizontal')
  
  %plot(hca,nanmean(vtrap_minus)*1e-3*[1 1],hca.YLim)
  errorbar(hca,nanmean(vtrap_minus)*1e-3,yloc_errobar,std(vtrap_minus,'omitnan')*1e-3,'horizontal')
  
  
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 

  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 0 % f instability
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [finst,params_inst] = mms_20170706_135303.get_f0(2);
    plot(hca,v_plot*v_scale,finst(v_plot*1e3,params_inst.n,params_inst.vd,params_inst.vt))
    str_info = {['T_{in}= [' sprintf('%g  ',params_inst.T) '] eV'];...
                ['n_{in}= [' sprintf('%g  ',params_inst.n*1e-6) '] cc'];...
                ['v_{d,in}= [' sprintf('%g  ',params_inst.vd*1e-3) '] km/s'];...
                };
    hleg = irf_legend(hca,str_info,[0.99,0.8]);
    hold(hca,'off')
end
hca.XLim = [-35 35];
h_lines = findobj(hca,'Type','line');
c_eval('h_lines(?).LineWidth = 1.5;',1:numel(h_lines))
hca.YLim = [0 2.6e-3];
set(hca,'ColorOrder',colors)
%irf_legend(hca,{'f_{FPI}';'f^{mod}';'f_0'},[0.02 0.99])
h_lines = findobj(hca,'type','line');
h_patches = findobj(hca,'type','patch');
legend([h_lines(end:-1:5); h_patches(end:-1:2); h_edi(1)],{'f_e^{FPI}';'f_e^{mod}';'f_0';'v_{ph}';'v_{trap}';'v_{EDI}'},...
    'location','northwest','box','off')
c_eval('h_edi(?).LineWidth = 1.0;',1:4)
hca.FontSize = 12;
end
if 0 % updated version for updated paper
  %%
colors = mms_colors('matlab');
v_fpi = mean(ef1D_bgremoved.depend{1},1);
ind0 = find(v_fpi==0);
doAverage = 1;
if doAverage
  f_fpi = mean(ef1D.data,1);
  f_fpi_bgrem = nanmean(ef1D_bgremoved.data,1);
  fpi_utc = ef1D.time([1 end]).utc;
  fpi_utc = sprintf('%s - %s',fpi_utc(1,12:23),fpi_utc(2,12:23));
else
  f_fpi = ef1D.data;
  f_fpi_bgrem = ef1D_bgremoved.data;
  fpi_utc = ef1D.time.utc;
end

f_fpi(:,ind0) = [];
f_fpi_bgrem(:,ind0) = [];
v_fpi(:,ind0) = [];

figure(113)
nrows = 1;
ncols = 1;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 0 % f fpi
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi,'-');
  hca.YLabel.String = {'f (s^1/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
  
  if doAverage
    str_lines = fpi_utc;
  else
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
  end
  irf_legend(hca,str_lines,[0.99 0.99])
  set(hca,'ColorOrder',zeros(10,3))
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % f fpi, rem bg
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi_bgrem,'-');
  hca.YLabel.String = {'f (s/m^4)'};
  hca.XLabel.String = {'v (10^3 km/s)'};  
  
  if doAverage
    str_lines = sprintf('FPI time interval:\n%s',fpi_utc);
  else
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
  end
  h_fpi_str = irf_legend(hca,str_lines,[0.99 0.99]);
  set(hca,'ColorOrder',zeros(10,3))
  h_fpi_str.HorizontalAlignment = 'left';
  h_fpi_str.Position(1) = 0.62;
  h_fpi_str.Position(2) = 0.95;
  h_fpi_str.BackgroundColor = [1 1 1];
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % <f> model
  hold(hca,'on')
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_vec*v_scale*1e-3,mean(Fabel_obs,1),'-','color',colors(2,:));
  hca.YLabel.String = {'f_e (s/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
      %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  %hca.XLim = 1.7*[-30 30];
  %hca.YLim = [0 3.5]*1e-3;
  hold(hca,'off')
end
if 1 % f0
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [f0,params_0] = mms_20170706_135303.get_f0(19);
    plot(hca,v_plot*v_scale,f0(v_plot*1e3,params_0.n,params_0.vd,params_0.vt),'color',colors(3,:))
    if 1
      str_info = {'f_0 input parameters:'
                ['T_{0}= [' sprintf('%g  ',params_0.T) '] eV'];...
                ['n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cm^{-3}'];...
                ['v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];...
                };
    elseif 0
      str_info = ['f_0 input parameters:\n ',...
                'T_{0}= [' sprintf('%g  ',params_0.T) '] eV\n ',...
                'n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cc\n ',...
                'v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];
    else
      str_info = sprintf('fpi time interval:\n%s',fpi_utc);
    end
    h_f0 = irf_legend(hca,str_info,[0.90,0.35]);
    c_eval('h_f0(?).HorizontalAlignment = ''left'';',1:4)
    c_eval('h_f0(?).Position(1) = 0.62;',1:4)
    %c_eval('h_f0(?).Position(2) = 0.7;',1:4)
  
    c_eval('h_f0(?).BackgroundColor = [1 1 1];',1:4)
    c_eval('h_f0(?).BackgroundColor = ''none'';',1:4)
    hold(hca,'off')
end
if 1 % v_edi
  hold(hca,'on')
  if 1 % plot dashed lines
    h_edi = plot(hca,[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim*0.20,'k--',-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--');
    %plot(hca,-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--')
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vph patch
  hold(hca,'on')
  all_vph = [mean(obs_velocity)-std(obs_velocity);...
             mean(obs_velocity);... 
             mean(obs_velocity)+std(obs_velocity)]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,all_vph*1e-3,hca.YLim,'k-.')
  else % plot patch                    
    hpatch = patch(hca,1e-3*[all_vph(1,1),all_vph(end,1) all_vph(end,1) all_vph(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.7 0.9]);
    hpatch.FaceAlpha = 0.5;
    hpatch.EdgeAlpha = 0.5;
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vtrap patch
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,plot_vtrap_minus*1e-3,hca.YLim,'k-.',plot_vtrap_plus*1e-3,hca.YLim,'k-.')
  else % plot patch
    hpatch = patch(hca,1e-3*[plot_vtrap_plus(1,1),plot_vtrap_plus(end,1) plot_vtrap_plus(end,1) plot_vtrap_plus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],colors(5,:)); % [0.8 0.8 0.8]
    hpatch.FaceAlpha = 0.3;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '-';

    hpatch = patch(hca,1e-3*[plot_vtrap_minus(1,1),plot_vtrap_minus(end,1) plot_vtrap_minus(end,1) plot_vtrap_minus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],colors(5,:)); % [0.8 0.8 0.8]
    hpatch.FaceAlpha = 0.3;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '-';
  end
  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
yloc_errobar = 2.3*1e-3;
if 0 % vph errorbar for std
  hold(hca,'on')  
  %plot(hca,mean(obs_velocity)*1e-3*[1 1],hca.YLim)
  errorbar(hca,mean(obs_velocity)*1e-3,yloc_errobar,std(obs_velocity)*1e-3,'horizontal')
  hold(hca,'off')
end
if 0 % vtrap errorbar for std
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  
  %plot(hca,nanmean(vtrap_plus)*1e-3*[1 1],hca.YLim,'-')
  errorbar(hca,nanmean(vtrap_plus)*1e-3,yloc_errobar,std(vtrap_plus,'omitnan')*1e-3,'horizontal')
  
  %plot(hca,nanmean(vtrap_minus)*1e-3*[1 1],hca.YLim)
  errorbar(hca,nanmean(vtrap_minus)*1e-3,yloc_errobar,std(vtrap_minus,'omitnan')*1e-3,'horizontal')
  
  
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 

  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 0 % f instability
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [finst,params_inst] = mms_20170706_135303.get_f0(2);
    plot(hca,v_plot*v_scale,finst(v_plot*1e3,params_inst.n,params_inst.vd,params_inst.vt))
    str_info = {['T_{in}= [' sprintf('%g  ',params_inst.T) '] eV'];...
                ['n_{in}= [' sprintf('%g  ',params_inst.n*1e-6) '] cc'];...
                ['v_{d,in}= [' sprintf('%g  ',params_inst.vd*1e-3) '] km/s'];...
                };
    hleg = irf_legend(hca,str_info,[0.99,0.8]);
    hold(hca,'off')
end
hca.XLim = [-35 35];
h_lines = findobj(hca,'Type','line');
c_eval('h_lines(?).LineWidth = 1.5;',1:numel(h_lines))
hca.YLim = [0 2.6e-3];
set(hca,'ColorOrder',colors)
%irf_legend(hca,{'f_{FPI}';'f^{mod}';'f_0'},[0.02 0.99])
h_lines = findobj(hca,'type','line');
h_patches = findobj(hca,'type','patch');
legend([h_lines(end:-1:5); h_patches(end:-1:2); h_edi(1)],{'<f_e^{FPI}>';'<f_e^{mod}>';'f_0';'v_{ph}';'v_{trap}';'v_{EDI}'},...
    'location','northwest','box','off')
c_eval('h_edi(?).LineWidth = 1.0;',1:4)
hca.FontSize = 12;

set(hca,'children',flipud(get(hca,'children')))        
end
if 0 % updated version for updated paper, different annotation
  %%
colors = mms_colors('matlab');
v_fpi = mean(ef1D_bgremoved.depend{1},1);
ind0 = find(v_fpi==0);
doAverage = 1;
if doAverage
  f_fpi = mean(ef1D.data,1);
  f_fpi_bgrem = nanmean(ef1D_bgremoved.data,1);
  fpi_utc = ef1D.time([1 end]) + [-0.015 0.015];
  fpi_utc = fpi_utc.utc;  
  fpi_utc = sprintf('%s - %s',fpi_utc(1,12:23),fpi_utc(2,12:23));
else
  f_fpi = ef1D.data;
  f_fpi_bgrem = ef1D_bgremoved.data;
  
  fpi_utc = ef1D.time.utc;
end

f_fpi(:,ind0) = [];
f_fpi_bgrem(:,ind0) = [];
v_fpi(:,ind0) = [];

figure(113)
nrows = 1;
ncols = 1;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 0 % f fpi
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi,'-');
  hca.YLabel.String = {'f (s^1/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
  
  if doAverage
    str_lines = fpi_utc;
  else
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
  end
  irf_legend(hca,str_lines,[0.99 0.99])
  set(hca,'ColorOrder',zeros(10,3))
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % f fpi, rem bg
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi_bgrem,'-');
  hca.YLabel.String = {'f (s/m^4)'};
  hca.XLabel.String = {'v (10^3 km/s)'};  
  
  if 0
    if doAverage
      str_lines = sprintf('FPI time interval:\n%s',fpi_utc);
    else
      str_lines = {...
        sprintf('-- fpi: %s',fpi_utc(1,12:23));...
        sprintf('-- fpi: %s',fpi_utc(2,12:23));...
        sprintf('-- fpi: %s',fpi_utc(3,12:23));...
        sprintf('-- fpi: %s',fpi_utc(4,12:23))};                  
    end
    h_fpi_str = irf_legend(hca,str_lines,[0.99 0.99]);
    set(hca,'ColorOrder',zeros(10,3))
    h_fpi_str.HorizontalAlignment = 'left';
    h_fpi_str.Position(1) = 0.62;
    h_fpi_str.Position(2) = 0.95;
    h_fpi_str.BackgroundColor = [1 1 1];
  end
  %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  hca.XLim = 1.7*[-30 30];
  hca.YLim = [0 3.5]*1e-3;
end
if 1 % <f> model
  hold(hca,'on')
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_vec*v_scale*1e-3,mean(Fabel_obs,1),'-','color',colors(2,:));
  hca.YLabel.String = {'f_e (s/m^4)'};
  hca.XLabel.String = {'v_{||} (10^3 km/s)'};  
      %irf_legend(hca,{sprintf('E_{low}=%g eV',unique(lowerelim.data))},[0.02 0.99])
  %hca.XLim = 1.7*[-30 30];
  %hca.YLim = [0 3.5]*1e-3;
  hold(hca,'off')
end
if 1 % f0
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [f0,params_0] = mms_20170706_135303.get_f0(19);
    plot(hca,v_plot*v_scale,f0(v_plot*1e3,params_0.n,params_0.vd,params_0.vt),'color',colors(3,:))
    if 0 % write f0 parmeters in figure
      if 1 % write f0 parmeters in figure
        str_info = {'f_0 input parameters:'
                  ['T_{0}= [' sprintf('%g  ',params_0.T) '] eV'];...
                  ['n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cm^{-3}'];...
                  ['v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];...
                  };
      elseif 0
        str_info = ['f_0 input parameters:\n ',...
                  'T_{0}= [' sprintf('%g  ',params_0.T) '] eV\n ',...
                  'n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cc\n ',...
                  'v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];
      else
        str_info = sprintf('fpi time interval:\n%s',fpi_utc);
      end
      h_f0 = irf_legend(hca,str_info,[0.90,0.35]);
      c_eval('h_f0(?).HorizontalAlignment = ''left'';',1:4)
      c_eval('h_f0(?).Position(1) = 0.62;',1:4)
      %c_eval('h_f0(?).Position(2) = 0.7;',1:4)

      c_eval('h_f0(?).BackgroundColor = [1 1 1];',1:4)
      c_eval('h_f0(?).BackgroundColor = ''none'';',1:4)
    end
    hold(hca,'off')
end
if 1 % v_edi
  hold(hca,'on')
  if 1 % plot dashed lines
    %h_edi = plot(hca,[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim*0.20,'k--',-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--');
    h_edi = plot(hca,[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--',-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--');
           %plot(hca,-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6,hca.YLim,'k--')
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vph patch
  hold(hca,'on')
  all_vph = [mean(obs_velocity)-std(obs_velocity);...
             mean(obs_velocity);... 
             mean(obs_velocity)+std(obs_velocity)]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,all_vph*1e-3,hca.YLim,'k-.')
  else % plot patch                    
    hpatch = patch(hca,1e-3*[all_vph(1,1),all_vph(end,1) all_vph(end,1) all_vph(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.7 0.9]);
    hpatch.FaceAlpha = 0.5;
    hpatch.EdgeAlpha = 0.5;
  end
  %hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vtrap patch
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,plot_vtrap_minus*1e-3,hca.YLim,'k-.',plot_vtrap_plus*1e-3,hca.YLim,'k-.')
  else % plot patch
    hpatch = patch(hca,1e-3*[plot_vtrap_plus(1,1),plot_vtrap_plus(end,1) plot_vtrap_plus(end,1) plot_vtrap_plus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],colors(5,:)); % [0.8 0.8 0.8]
    hpatch.FaceAlpha = 0.3;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '-';

    hpatch = patch(hca,1e-3*[plot_vtrap_minus(1,1),plot_vtrap_minus(end,1) plot_vtrap_minus(end,1) plot_vtrap_minus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],colors(5,:)); % [0.8 0.8 0.8]
    hpatch.FaceAlpha = 0.3;
    hpatch.EdgeAlpha = 0.3;
    hpatch.LineStyle = '-';
  end
  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
yloc_errobar = 2.3*1e-3;
if 0 % vph errorbar for std
  hold(hca,'on')  
  %plot(hca,mean(obs_velocity)*1e-3*[1 1],hca.YLim)
  errorbar(hca,mean(obs_velocity)*1e-3,yloc_errobar,std(obs_velocity)*1e-3,'horizontal')
  hold(hca,'off')
end
if 0 % vtrap errorbar for std
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  
  %plot(hca,nanmean(vtrap_plus)*1e-3*[1 1],hca.YLim,'-')
  errorbar(hca,nanmean(vtrap_plus)*1e-3,yloc_errobar,std(vtrap_plus,'omitnan')*1e-3,'horizontal')
  
  %plot(hca,nanmean(vtrap_minus)*1e-3*[1 1],hca.YLim)
  errorbar(hca,nanmean(vtrap_minus)*1e-3,yloc_errobar,std(vtrap_minus,'omitnan')*1e-3,'horizontal')
  
  
  plot_vtrap_plus = [nanmean(vtrap_plus)-std(vtrap_plus,'omitnan');...
                     nanmean(vtrap_plus);...
                     nanmean(vtrap_plus)+std(vtrap_plus,'omitnan')]*[1 1]; 
  plot_vtrap_minus = [nanmean(vtrap_minus)-std(vtrap_minus,'omitnan');...
                      nanmean(vtrap_minus);...
                      nanmean(vtrap_minus)+std(vtrap_minus,'omitnan')]*[1 1]; 

  %hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';

  %hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  %hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 0 % f instability
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [finst,params_inst] = mms_20170706_135303.get_f0(2);
    plot(hca,v_plot*v_scale,finst(v_plot*1e3,params_inst.n,params_inst.vd,params_inst.vt))
    str_info = {['T_{in}= [' sprintf('%g  ',params_inst.T) '] eV'];...
                ['n_{in}= [' sprintf('%g  ',params_inst.n*1e-6) '] cc'];...
                ['v_{d,in}= [' sprintf('%g  ',params_inst.vd*1e-3) '] km/s'];...
                };
    hleg = irf_legend(hca,str_info,[0.99,0.8]);
    hold(hca,'off')
end
hca.XLim = [-35 35];
h_lines = findobj(hca,'Type','line');
c_eval('h_lines(?).LineWidth = 1.5;',1:numel(h_lines))
hca.YLim = [0 2.6e-3];
set(hca,'ColorOrder',colors)
%irf_legend(hca,{'f_{FPI}';'f^{mod}';'f_0'},[0.02 0.99])
h_lines = findobj(hca,'type','line');
h_patches = findobj(hca,'type','patch');
legend([h_lines(end:-1:5); h_patches(end:-1:2); h_edi(1)],{'<f_e^{FPI}>';'<f_e^{mod}>';'f_0';'v_{ph}';'v_{trap}';'v_{EDI}'},...
    'location','northeast','box','off')
c_eval('h_edi(?).LineWidth = 1.0;',1:4)
hca.FontSize = 12;

set(hca,'children',flipud(get(hca,'children')))  

hca.Title.String = ['FPI time interval:' fpi_utc];
hca.Title.FontWeight = 'normal';
hca.XLim = [-25 25];

end