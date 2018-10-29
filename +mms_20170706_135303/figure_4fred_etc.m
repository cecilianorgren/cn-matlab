units = irf_units;
% Load datastore
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   

% MMS id and time interval
mms_id = 1;
tint = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.620Z')+3;
c_eval('eDist = ePDist?;',mms_id)

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

% Remove background electrons
nPhoto = 00;
nSecond = 20;
if 0
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint, ...
             'Nphotoe_art', nPhoto, 'nSecondary', nSecond, 'ZeroNaN', 0);
else
[eDist_bgremoved, eDist_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tint,...
             'ZeroNaN', 0);  
end
%% Make reduced distribution
strTint = [irf_time(tint(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vint = [-Inf Inf];
vg = (-70:2:70)*1e3;
c_eval('eDist = ePDist?.tlim(tint);',mms_id)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 00;
tic; ef1D = eDist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B
tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc 
tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc
tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc

tic; ef1D_bgremoved = eDist_bgremoved.reduce('1D',ePara,'vint',vint,'scpot',scpot); toc % reduced distribution along B
tic; ef2D_parperp1_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc 
tic; ef2D_parperp2_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc
tic; ef2D_perp1perp2_bgremoved = eDist_bgremoved.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart'); toc

%%
for it = 1%1:4
  %%
  figure
  %it = 3;
  nrows = 2;
  ncols = 4;
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
  
  h(1).Title.String = {ef2D_parperp1(it).time.utc,'no background removed'};
  h(5).Title.String = {'background removed'};
  for ipanel = [1 5]
    h(ipanel).CLim = [-13 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg([1 end])*1e-3*0.5;
    h(ipanel).YLim = [0 4e-3];
  end  
  for ipanel = [2:4 6:8]
    h(ipanel).CLim = [-13 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg([1 end])*1e-3*1;
    h(ipanel).YLim = vg([1 end])*1e-3*1;
  end
  %cn.print(sprintf('f1D_2D_bgrem_it=%g_nPhoto=%g_nSecond=%g',it,nPhoto,nSecond))
end

% Run instability anaysis
[f_ins,params_ins] = mms_20170706_135303.get_f0(2);
[f0,params_0] = mms_20170706_135303.get_f0(1);

%% 4 fred, vph/trapping range, EDI range, average flux from model, and f0 of model

if 0
  %%
figure
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

if 1 % f   
  hca = h(isub); isub = isub + 1;
  v_scale = 1e-3;
  %hlines = plot(hca,v_fpi*v_scale,f_fpi,'-',v_fpi*v_scale,f_fpi_bgrem,'--');
  hlines = plot(hca,v_fpi*v_scale,f_fpi_bgrem,'-');
  hca.YLabel.String = {'f (s^1/m^4)'};
  hca.XLabel.String = {'v (10^3 km/s)'};  
  
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
if 1 % vph 
  hold(hca,'on')
  all_vph = [mean(obs_velocity)-std(obs_velocity);...
             mean(obs_velocity);...
             mean(obs_velocity)+std(obs_velocity)]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,all_vph*1e-3,hca.YLim,'k-.')
  else % plot patch                    
    hpatch = patch(hca,1e-3*[all_vph(1,1),all_vph(end,1) all_vph(end,1) all_vph(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.8 0.8]);
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeAlpha = 0.3;
  end
  hleg = irf_legend(hca,'vph',[0.5*21/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % vtrap
  hold(hca,'on')
  obs_vtrap = obs_vtrap_all(:,mms_id);
  vtrap_plus = obs_velocity + obs_vtrap*1e-3;
  vtrap_minus = obs_velocity - obs_vtrap*1e-3;
  plot_vtrap_plus = [mean(vtrap_plus)-std(vtrap_plus);...
                     mean(vtrap_plus);...
                     mean(vtrap_plus)+std(vtrap_plus)]*[1 1]; 
  plot_vtrap_minus = [mean(vtrap_minus)-std(vtrap_minus);...
                      mean(vtrap_minus);...
                      mean(vtrap_minus)+std(vtrap_minus)]*[1 1]; 
  if 0 % plot dashed lines
    plot(hca,plot_vtrap_minus*1e-3,hca.YLim,'k-.',plot_vtrap_plus*1e-3,hca.YLim,'k-.')
  else % plot patch
    hpatch = patch(hca,1e-3*[plot_vtrap_plus(1,1),plot_vtrap_plus(end,1) plot_vtrap_plus(end,1) plot_vtrap_plus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.8 0.8]);
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeAlpha = 0.3;

    hpatch = patch(hca,1e-3*[plot_vtrap_minus(1,1),plot_vtrap_minus(end,1) plot_vtrap_minus(end,1) plot_vtrap_minus(1,1)],...
      [0 0 hca.YLim(2) hca.YLim(2)],[0.8 0.8 0.8]);
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeAlpha = 0.3;
  end
  hleg = irf_legend(hca,'vph-vtrap',[0.5*12/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  hleg.HorizontalAlignment = 'center';

  hleg = irf_legend(hca,'vph+vtrap',[0.5*30/hca.XLim(2) 0.95],[0 0 0]);
  %hleg.BackgroundColor = [1 1 1];
  hleg.HorizontalAlignment = 'center';
  hold(hca,'off')
end
if 1 % f0
    hold(hca,'on')
    v_plot = vg(1):500:vg(end);
    [f0,params_0] = mms_20170706_135303.get_f0(1);
    plot(hca,v_plot*v_scale,f0(v_plot*1e3,params_0.n,params_0.vd,params_0.vt))
    str_info = {['T_{0}= [' sprintf('%g  ',params_0.T) '] eV'];...
                ['n_{0}= [' sprintf('%g  ',params_0.n*1e-6) '] cc'];...
                ['v_{d,0}= [' sprintf('%g  ',params_0.vd*1e-3) '] km/s'];...
                };
    hleg = irf_legend(hca,str_info,[0.99,0.1]);
    hold(hca,'off')
end
if 1 % f instability
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
end