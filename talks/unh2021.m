%% Movie of streaming instability and thermalization
localuser = datastore('local','user');
smilei_path = datastore('local','smilei_data_dir');
directory = [smilei_path '/streaming_instabilities/cn_beam_1d/beam1c/'];
namelist = [directory 'tst1d_02_two_str_instability.py'];

%directory = ['/Users/cno062/Discs/betzy/Smilei/cn_beam_1d/beam_buneman/'];
%namelist = [directory 'buneman.py'];

filepath = [directory 'Fields0.h5'];
particlebinningpath = directory;

%sm = SMILEI(filepath,namelist,particlebinningpath);

info = h5info(filepath);
nGroups = numel(info.Groups.Groups);
datasets = {info.Groups.Groups(1).Datasets.Name};
nDatasets = numel(datasets);
for iDataset = 1:nDatasets
  clear(datasets{iDataset});
end

for iGroup = 1:nGroups
  group = info.Groups.Groups(iGroup).Name;
  for iDataset = 1:nDatasets
    dataset = h5read(filepath,[group filesep datasets{iDataset}]);
    eval([datasets{iDataset} '(iGroup,:) = torow(dataset);' ])    
  end
end
x0 = h5readatt(filepath,[group filesep datasets{iDataset}],'gridGlobalOffset');
dx = h5readatt(filepath,[group filesep datasets{iDataset}],'gridSpacing');
nx = info.Groups.Groups(1).Datasets(1).Dataspace.Size;
x = linspace(x0,nx*dx,nx);

particlebinningpath = directory;
filename = sprintf('ParticleBinning%g.h5',0);
info_diag = h5info([particlebinningpath filesep filename]);
datasets = {info_diag.Datasets.Name};
nDatasets = numel(datasets);

%% three panels, Ex, rho, f(x,vx)
nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols);
h(1).Position = [0.15 0.7 0.7 0.2];
h(2).Position = [0.15 0.5 0.7 0.2];
h(3).Position = [0.15 0.15 0.7 0.35];

clim = [0 3e-4];
doLog = 0;

for it = 1:nGroups
  isub = 1;

  if 1 % Ex
    hca = h(isub); isub = isub + 1;
    plot(hca,x,Ex(it,:))
    hca.YLim = 0.1*[-1 1];
    hca.XLim = x([1 end]);
  end
  
  if 1 % Rho
    hca = h(isub); isub = isub + 1;
    plot(hca,x,Rho(it,:))
    hca.YLim = 1*[-1 1];
    hca.XLim = x([1 end]);
  end 

  if 1 % f(x,vx)
    hca = h(isub); isub = isub + 1;
    dataset = info_diag.Datasets(it).Name; 
    data = h5read([particlebinningpath filename],[filesep dataset]);
    if it == 1
      data0 = data;
    end 
    name = h5readatt([particlebinningpath filename],'/','name');
    dep0 = h5readatt([particlebinningpath filename],'/','axis0');
    dep1 = h5readatt([particlebinningpath filename],'/','axis1');  
    [dep0_str,dep0_val] = make_dep(dep0);
    [dep1_str,dep1_val] = make_dep(dep1); 
    if doLog
      imagesc(hca,dep0_val,dep1_val,smooth2(log10(squeeze(data)),2))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = ['log10 ' name];
    else
      imagesc(hca,dep0_val,dep1_val,smooth2(squeeze(data),2))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = name;
    end
    %hca.CLim = clim;
    hca.XLabel.String = dep0_str;
  %    hca.YLabel.String = dep2_str;
    hca.Title.String = dataset;
    hcb.YLabel.Interpreter = 'none';
    hca.YDir = 'normal';
    colormap(hca,pic_colors('candy4'))
    hca.XLim = x([1 end]);
  end
  if 0 % <f(x)>
    hca = h(isub); isub = isub + 1;
    if doLog
      plot(dep1_val,mean(log10(squeeze(data)),2),dep1_val,mean(log10(squeeze(data0)),2))
      hca.YLabel.String = ['log10 ' name];    
    else
      plot(dep1_val,mean(squeeze(data),2),dep1_val,mean(squeeze(data0),2))
      hca.YLabel.String = name;
    end
    
    hca.XLabel.String = dep0_str;
    hca.Title.String = dataset;
    hca.YLabel.Interpreter = 'none';   
  end
  
  %compact_panels
  irf_plot_axis_align
  drawnow
  pause(0.1)
end

%% 2 panels, 1: (Ex, rho, f(x,vx)), 2: <f(x,vx)>
% Video stuff
doVideo = 1;
fileName = [printpath 'test1'];

fig = figure(123);
fig.Color = [1 1 1];

nrows = 1;
ncols = 2;
h = setup_subplots(nrows,ncols);
h(1).Position = [0.10 0.15 0.50 0.75];
h(2).Position = [0.70 0.15 0.25 0.75];

fontsize = 13;
clim = [0 1e-4];
doLog = 0;
vlim = 0.3*[-1 1];


% Allow for adjustments
disp('Adjust figure size, then hit any key to continue.')
pause

% Setup for gif an video
if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 5;
  open(vidObj);        
end

 % electron hole ends up at the border of the axes, shifts it so that it ends up in the middle
doCircshift = 1; 
nCircshift = 200;
nCircshiftFields = fix(size(Ex,2)/2);

nSmooth = 1;

for it = 1:nGroups
  isub = 1;
 

  if 1 % f(x,vx)
    hca = h(isub); isub = isub + 1;
    dataset = info_diag.Datasets(it).Name; 
    data = h5read([particlebinningpath filename],[filesep dataset]);
    if it == 1
      data0 = data;
    end 
    
    if doCircshift
      data = circshift(data,nCircshift,2);
    end
    
    name = h5readatt([particlebinningpath filename],'/','name');
    dep0 = h5readatt([particlebinningpath filename],'/','axis0');
    dep1 = h5readatt([particlebinningpath filename],'/','axis1');  
    [dep0_str,dep0_val] = make_dep(dep0);
    [dep1_str,dep1_val] = make_dep(dep1); 
    if doLog
      imagesc(hca,dep0_val,dep1_val,smooth2(log10(squeeze(data)),nSmooth))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = ['log10 ' name];
      hca.CLim = [1e-3 log10(clim(2))];
    else
      imagesc(hca,dep0_val,dep1_val,smooth2(squeeze(data),nSmooth))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = name;
      hca.CLim = clim;
    end
    
    hca.YLim = vlim;
    hca.YLabel.String = 'v_x';
    hca.XLabel.String = dep0_str;
  %    hca.YLabel.String = dep2_str;
    %hca.Title.String = dataset;
    %hcb.YLabel.Interpreter = 'none';
    hcb.YLabel.String = 'f(x,v_x)';
    hcb.YLabel.FontSize = fontsize;
    hca.YDir = 'normal';
    colormap(hca,pic_colors('candy4'))
    hca.XLim = x([1 end]);
    
    hold(hca,'on')
    if 1 % Ex      
      if doCircshift
        plot(hca,x,3*circshift(Ex(it,:),nCircshiftFields,2),'k','linewidth',2)
      else
        plot(hca,x,2*Ex(it,:),'k','linewidth',1)
      end
      %hca.YLim = 0.1*[-1 1];
      hca.XLim = x([1 end]);
    end
    if 1 % Rho           
      if doCircshift
        plot(hca,x,0.4*circshift(Rho(it,:),nCircshiftFields,2),'k-.','linewidth',2)
      else
        plot(hca,x,0.4*Rho(it,:),'k','linewidth',1)
      end
      %hca.YLim = 1*[-1 1];
      hca.XLim = x([1 end]);
    end 
    hold(hca,'off')
    
    
    %hl = findobj(h(1),'type','line');
    %legend(hl(2:-1:1),{'E_{||}','\rho'},'box','off','location','northwest')
    
    
  end
  if 1 % <f(x)>
    hca = h(isub); isub = isub + 1;
    if doLog
      plot(dep1_val,mean(log10(squeeze(data)),2),dep1_val,mean(log10(squeeze(data0)),2),'k-.','linewidth',1)
      hca.YLabel.String = ['log10 ' name];    
    else
      plot(dep1_val,mean(squeeze(data),2),dep1_val,mean(squeeze(data0),2),'k-.','linewidth',1)
      hca.YLabel.String = name;
    end
    
    hca.XLabel.String = 'v_x';    
    hca.YLabel.String = '<f(x,v_x)>';
    %hca.Title.String = dataset;
    %hca.YLabel.Interpreter = 'none';   
    hca.XLim = vlim;
    
  end
  
  %compact_panels
  irf_plot_axis_align
  for ip = 1:ncols
    h(ip).FontSize = fontsize;
  end
  
  % Overall formatting
  h(1).Title.String = 'Electron phase space';
  h(2).Title.String = 'Averaged phase space';
  hcb.YTick = [];
  hcb.YLabel.String = [];
  %h(1).YTick = [];
  h(2).YTick = [];
  h(1).YTick = -0.2:0.1:0.2;
  h(2).XTick = -0.2:0.1:0.2;
  %h(1).XTick = [];
  %h(2).XTick = [];  
  hcb.YLabel.String = 'f(x,v_x)';
  h(2).YLabel.String = 'f(v_x)';
  h(1).Position(3) = 0.50;
  h(2).Position(1) = 0.72;
  h(1).YLabel.String = 'v_x/c';
  h(1).XLabel.String = 'x/d_e';
  h(2).XLabel.String = 'v_x/c';  
  
  drawnow
  pause(0.1)
  if doVideo
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
end

% Write gif and video
if doVideo
  close(vidObj)
end


%% Overview pictures, no02m
%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapth = pic_colors('thermal');
cmapca = pic_colors('candy4');

twpe = 20000;
xlim = no02m.xi([1 end]);
zlim = no02m.zi([1 end]);
xlim = [80 105];
zlim = [-7 7];

varstrs = {'vix','vex','ni','Epar'}';
clims = {[-1 1],[-10 10],[0 1],[-1 1]};
cmaps = {cmapbr,cmapbr,cmapth,cmapbr};


varstrs = {'vix','vex','Epar'}';
clims = {[-1 1]*0.99,[-10 10]*0.99,0.8*[-1 1]};
cmaps = {cmapbr,cmapbr,cmapbr};
cbarlabels = {'v_{ix}','v_{ex}','E_{||}'};

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'smooth',2,'cbarlabels',cbarlabels);

%% MMS ESW overview
% mms_20170706_135303.load_data
% mms_20170706_135303.prepare_data_single_sc

% first part of mms_20170706_135303.paper_fig1.m

%% Plot, zoom only, edi fpi comparison
ic = 1;
npanels = 3;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);

% common time shifts
dt = [0.0000  -0.0012  -0.0009  -0.0012];
dt0 = 0.0008;
dt = dt + dt0;

doDt = 0;
if not(doDt), dt = [0,0,0,0]; end
doFour = 0;
doOne = 1;
  
if 0 % E par, 4 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if doFour % E par, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},...
                 [0.01 1.05],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if doOne % E par, single spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if doFour % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if doOne % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('hh = irf_plot(hca,{phi?},''comp'',''dt'',dt);',ic)
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
end
if doFour % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_divE/units_scaling,...
                     },'comp','dt',[dt 0]); 
  hh(1).Children(1).LineWidth = 2;
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
    set(hca,'ColorOrder',mms_colors('1234b'))
    irf_legend(hca,{'mms1';'mms2';'mms3';'mms4';'4 sc'},[1.02 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'4 sc'},[1.02 0.9],'fontsize',12);
  end
end
if doOne % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  if doFour
    c_eval('hh = irf_plot(hca,{dn_E?par/v_for_density_scaling*1e-6/units_scaling},''comp'',''dt'',[dt 0]);',ic) 
    hh(1).Children(1).LineWidth = 2;
  else
    c_eval('irf_plot(hca,{dn_E?par/v_for_density_scaling*1e-6/units_scaling},''comp'',''dt'',[dt 0]);',ic)     
  end
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if doFour % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';  
  end
end

%hca = irf_panel('edi flux'); hca.YLim = [0 7.999];

irf_zoom(h,'x',phi1.time)
irf_zoom(h,'y')
irf_plot_axis_align


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

hca = irf_panel('E par'); hca.YLim = [-60 60];
hca = irf_panel('phi'); hca.YLim = [0.001 399];
hca = irf_panel('density perturbation'); hca.YLim = [-1 1.8];


doDoubleAxis = 0; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
  set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
  set(ax2,'YAxisLocation','right');
  set(ax2,'Color','none','box','off'); % color of axis      
  %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
  %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
  irf_timeaxis(ax2,'nolabels')
  ax2.XLabel.String = [];
  ax2.YLabel.String = {'\delta n/n'};
  ax2.YLabel.Interpreter = 'tex';    
  ax2.YTick = hca.YTick*1e-3/n0;  
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align(h)
%ax.YLabel.Position(1) = 1.07;

%set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required

%% Plot, zoom only, Epar edi comparison
ic = 4;
npanels = 2;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);

% common time shifts
dt = [0.0000  -0.0012  -0.0009  -0.0012];
dt0 = 0.0008;
dt = dt + dt0;

doDt = 0;
if not(doDt), dt = [0,0,0,0]; end
doFour = 0;
doOne = 1;
  
if 0 % E par, 4 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0*doFour % E par, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},...
                 [0.01 1.05],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if doOne % E par, single spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 0*doFour % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 0*doOne % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('hh = irf_plot(hca,{phi?},''comp'',''dt'',dt);',ic)
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
end
if 0*doFour % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_divE/units_scaling,...
                     },'comp','dt',[dt 0]); 
  hh(1).Children(1).LineWidth = 2;
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
    set(hca,'ColorOrder',mms_colors('1234b'))
    irf_legend(hca,{'mms1';'mms2';'mms3';'mms4';'4 sc'},[1.02 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'4 sc'},[1.02 0.9],'fontsize',12);
  end
end
if 0*doOne % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  if doFour
    c_eval('hh = irf_plot(hca,{dn_E?par/v_for_density_scaling*1e-6/units_scaling},''comp'',''dt'',[dt 0]);',ic) 
    hh(1).Children(1).LineWidth = 2;
  else
    c_eval('irf_plot(hca,{dn_E?par/v_for_density_scaling*1e-6/units_scaling},''comp'',''dt'',[dt 0]);',ic)     
  end
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if doFour % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';  
  end
end
if 1 % edi flux 0 180 1 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(180)*1e-6},''comp'',''dt'',dt);',ic)
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}m^{-2} sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
end

%hca = irf_panel('edi flux'); hca.YLim = [0 7.999];

irf_zoom(h,'x',phi1.time)
irf_zoom(h,'y')
hca = irf_panel('edi flux'); hca.YLim = [0 6.99];
irf_plot_axis_align
%%

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

hca = irf_panel('E par'); hca.YLim = [-60 60];
hca = irf_panel('phi'); hca.YLim = [0.001 399];
hca = irf_panel('density perturbation'); hca.YLim = [-1 1.8];


doDoubleAxis = 0; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
  set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
  set(ax2,'YAxisLocation','right');
  set(ax2,'Color','none','box','off'); % color of axis      
  %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
  %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
  irf_timeaxis(ax2,'nolabels')
  ax2.XLabel.String = [];
  ax2.YLabel.String = {'\delta n/n'};
  ax2.YLabel.Interpreter = 'tex';    
  ax2.YTick = hca.YTick*1e-3/n0;  
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align(h)
%ax.YLabel.Position(1) = 1.07;

%set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required

%% Anticorrelation of phi and jedi

h = setup_subplots(1,1);
hca = h(1);
xx = phi1.resample(ePitch1_flux_edi).tlim(phi1.time).data;
yy = ePitch1_flux_edi.palim(palim).tlim(phi1.time).data*1e-6;
plot(xx,yy,'.k','markersize',15)
hold on
xvec = linspace(0,400,100);
%hfit = plot(xvec,(1/0.007)./xvec,'linewidth',1);
hfit = plot(xvec,(140)./xvec,'linewidth',1);
%scatter(xx,yy,'.')
hca.YLim = [0 4];
hca.XLim = [0 350];
hca.XLabel.String = '\phi (V)';
hca.YLabel.String = 'j_e^{EDI} (10^6 s^{-1}cm^{-2}sr^{-1})';
hca.FontSize = 12;

%% Anticorrelation of phi and jedi

h = setup_subplots(1,1);
hca = h(1);
xx = phi1.resample(ePitch1_flux_edi).tlim(phi1.time).data;
yy = ePitch1_flux_edi.palim(palim).tlim(phi1.time).data*1e-6;
yy = 1./yy;
yy(yy==Inf) = NaN;
plot(xx,yy,'.k','markersize',15)
hold on
xvec = linspace(0,400,100);
%hfit = plot(xvec,50./xvec,'linewidth',1);
%scatter(xx,yy,'.')
hold off
hca.YLim = [0 4];
hca.XLim = [0 350];
hca.XLabel.String = '\phi (V)';
hca.YLabel.String = 'j_e^{EDI} (10^6 s^{-1}cm^{-2}sr^{-1})';
hca.FontSize = 12;

%% More EDI plots

%% Plot, Epar, phi, dn, jedi
ic = 1;
npanels = 4;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);

% load eh data
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_vtrap = sqrt(2*units.e*obs_potential/units.me)*1e-3; % km/s
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
c_eval('obs_vtrap? = irf.ts_scalar(obs_t0_epoch_mms?,obs_vtrap(:,?));')
c_eval('obs_vph?.data(isnan(obs_potential(:,?))) = NaN;')


% common time shifts
dt = [0.0000  -0.0012  -0.0009  -0.0012];
dt0 = 0.0008;
dt = dt + dt0;
  
if 1 % E par, 4 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0%1 % E par, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 0 % E perp, abs, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1perp.abs,gseE2perp.abs,gseE3perp.abs,gseE4perp.abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 0%1 % E perp, abs, filt, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).abs,gseE2perp.filt(flow,0,[],5).abs,gseE3perp.filt(flow,0,[],5).abs,gseE4perp.filt(flow,0,[],5).abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  %irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f > %g Hz',flow)},[0.05 0.98],'fontsize',12);
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 0 % E perp, z, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp z dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).z,gseE2perp.filt(flow,0,[],5).z,gseE3perp.filt(flow,0,[],5).z,gseE4perp.filt(flow,0,[],5).z},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f>%g Hz',flow)},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,z}','(mV/m)'};  
end
if 0 % E perp, y, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp y dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).y,gseE2perp.filt(flow,0,[],5).y,gseE3perp.filt(flow,0,[],5).y,gseE4perp.filt(flow,0,[],5).y},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f>%g Hz',flow)},[0.98 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};  
end
if 0 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{obs_phi1,obs_phi2,obs_phi3,obs_phi4},'comp');
  %c_eval('hh.Children(?).Marker = ''.'';',1:4)
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 1 % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 1 % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
%   hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
%                      dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
%                      dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
%                      dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
%                      dn_divE/units_scaling,...
%                      },'comp','dt',[dt 0]); 
  hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
                     },'comp','dt',[dt]); 
  %hh(1).Children(1).LineWidth = 2;
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 0% 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
    set(hca,'ColorOrder',mms_colors('1234b'))
    irf_legend(hca,{'mms1';'mms2';'mms3';'mms4';'4 sc'},[1.02 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'4 sc'},[1.02 0.9],'fontsize',12);
  end
end
if 0 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fred vph vtrap');
  
  [hsurf,hcbar] = irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  hold(hca,'on')
  c_eval('vmin? = obs_vph? - obs_vtrap?;',1:4)
  c_eval('vmax? = obs_vph? + obs_vtrap?;',1:4)
  
  %set(hca,'ColorOrder',mms_colors('111223344'))
  set(hca,'ColorOrder',mms_colors('111223344'))
  vscale = 1e-3;
  hlines = irf_plot(hca,{obs_vph1*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  htrap = hlines.Children(1:end-2);
  hvph = hlines.Children(end-1);
  c_eval('htrap(?).Marker = ''.'';',1:numel(htrap))
  hvph.LineStyle = '--';
  
  %irf_patch(hca,{vmin,vmax})
  %hca.YLim = sort(real([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]));
  hold(hca,'off')
  hca.YLabel.String = {'v_{||}','(10^3 km/s)'};  
  if 0 % label vtrap vph vtrap
    set(hca,'ColorOrder',mms_colors('122'))
    irf_legend(hca,{'v_{ph}'},[0.4 0.7],'fontsize',12);
    irf_legend(hca,{'v_{trap}'},[0.4 0.99],'fontsize',12);
    irf_legend(hca,{'v_{trap}'},[0.4 0.3],'fontsize',12);
  end
  if 1 % label -- vph
    set(hca,'ColorOrder',mms_colors('1'))
    irf_legend(hca,{'-- v_{ph}'},[0.01 0.98],'fontsize',12);    
  end
  hsurf = findobj(hca.Children,'Type','Surface');
  hsurf.FaceAlpha = 1;
  hline = findobj(hca.Children,'Type','Line');
  c_eval('hline(?).LineWidth = 1.5;',1:numel(hline))
end
if 0 % edi flux 0 180 1 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{flux0_mms?,flux180_mms?},''comp'',''dt'',dt);',ic)
  hca.YLabel.String = {'flux 180^o','10^6 s^{-1}m^{-2})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
end
if 0 % edi flux 0 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux par');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [0 12];  
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6*0.4,ePitch3_flux_edi.palim(palim)*1e-6*0.5,ePitch4_flux_edi.palim(palim)*1e-6*0.5},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','parallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [0 11.25]^o'},[0.05 0.99],'fontsize',12);
end
if 1 % edi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [168 180];  
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6*0.4,ePitch3_flux_edi.palim(palim)*1e-6*0.5,ePitch4_flux_edi.palim(palim)*1e-6*0.5},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'\theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end

if 0 % edi flux 180 4sc, averaged to fpi and stepped
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux apar step');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi_apar_step*1e-6,ePitch2_flux_edi_apar_step*1e-6,ePitch3_flux_edi_apar_step*1e-6,ePitch4_flux_edi_apar_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [168.25, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'averaged to FPI timeline'},[0.98 0.98],'fontsize',12,'color',[0 0 0]);
end

if 0 % fpi flux 0 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi_par500_step*1e-6,ePitch2_flux_fpi_par500_step*1e-6,ePitch3_flux_fpi_par500_step*1e-6,ePitch4_flux_fpi_par500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','parallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end
if 0 % fpi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi_apar500_step*1e-6,ePitch2_flux_fpi_apar500_step*1e-6,ePitch3_flux_fpi_apar500_step*1e-6,ePitch4_flux_fpi_apar500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [168.25, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
  
end
if 0 % fpi flux 180 4sc, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar 2');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi2_apar500_step*1e-6,ePitch2_flux_fpi2_apar500_step*1e-6,ePitch3_flux_fpi2_apar500_step*1e-6,ePitch4_flux_fpi2_apar500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [157.50, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
end
if 0 % fpi flux 180 4sc, 11.25, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar 12 comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  hlines1 = irf_plot(hca,{ePitch1_flux_fpi_apar500_step*1e-6,ePitch2_flux_fpi_apar500_step*1e-6,ePitch3_flux_fpi_apar500_step*1e-6,ePitch4_flux_fpi_apar500_step*1e-6},'comp');
  hcl = hca.Children;  
  c_eval('hcl(?).LineWidth = 1;',1:4)
  hold(hca,'on')
  hlines2 = irf_plot(hca,{ePitch1_flux_fpi2_apar500_step*1e-6,ePitch2_flux_fpi2_apar500_step*1e-6,ePitch3_flux_fpi2_apar500_step*1e-6,ePitch4_flux_fpi2_apar500_step*1e-6},'comp');
  hcl = hca.Children;  
  c_eval('hcl(?).LineStyle = ''--'';',5:8)
  c_eval('hcl(?).LineWidth = 0.5;',5:8)
  c_eval('hcl(?).LineWidth = 1;',1:4)
  hold(hca,'off')
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'dashed - \theta = [168.75, 180]^o','    solid - \theta = [157.50, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
end
if 0 % jedi-jfpe fpi flux 180 4sc, 11.25, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi edi diff');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('dj? = ePitch?_flux_edi_apar_step+-ePitch?_flux_fpi_apar500_step;',1:4)  
  hlines1 = irf_plot(hca,{dj1*1e-6,dj2*1e-6,dj3*1e-6,dj4*1e-6},'comp');
  hcl = hca.Children;  
  %c_eval('hcl(?).LineWidth = 1;',1:4)  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}-j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'EDI - \theta = [168.75, 180]^o','    FPI - \theta = [168.75, 180]^o'},[0.05 0.02],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.98 0.98],'fontsize',12,'color',[0 0 0]);
end
if 0 % jedi-jfpe fpi flux 180 4sc, 11.25, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi edi diff larger range');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('dj? = ePitch?_flux_edi_apar_step+-ePitch?_flux_fpi2_apar500_step;',1:4)  
  hlines1 = irf_plot(hca,{dj1*1e-6,dj2*1e-6,dj3*1e-6,dj4*1e-6},'comp');
  hcl = hca.Children;  
  %c_eval('hcl(?).LineWidth = 1;',1:4)  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}-j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'EDI - \theta = [168.75, 180]^o','    FPI - \theta = [157.50, 180]^o'},[0.05 0.02],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.99 0.98],'fontsize',12,'color',[0 0 0]);
end


if 0 % edi, phi comparison for one spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi phi comp');
  set(hca,'ColorOrder',mms_colors('a1'))
  c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6},''dt'',dt(1),''color'',[0 0 0]);',ic)
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6});
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
  %ax2 = axis('position',hca.Position);
  hca.YLim = [0 3.9999];
  yyaxis('right')
  ax = gca;
  irf_plot(ax,phi1,'dt',dt(1))
  ax.YLabel.String = {'\phi (V)'};
  ax.YLabel.Interpreter = 'tex';
  %ax.YLabel.Position(1) = 1.07;
  irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'mms1'},[0.98 0.99],'fontsize',12,'color',[0 0 0]);
end
%c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
%hcb_fred.Position(2) = hcb_fred.Position(2)+0.02;
%c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',isub_long)


%irf_zoom(h(isub_long),'x',tint)
%irf_zoom(h(isub_short),'x',tint_zoom)
irf_zoom(h,'x',phi1.time)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(isub_long(end)).XLabel.String = [];

%[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(isub_long(end)),h(isub_short(1))); 
if 0 % two black lines
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom(1),[0 0 0]);',isub_long)
  c_eval('h_mark2(?) = irf_pl_mark(h(?),tint_zoom(2),[0 0 0]);',isub_long)
elseif 0 % colored region
  mark_color = [1 0.3 0];
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_long)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long(1:2))
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

%colormap(cn.cmap('blue_white'))
%colormap(cn.cmap('white_blue'))
%hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
%hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
%hca = irf_panel('Vi'); hca.YLim = [-799 399];
hca = irf_panel('edi flux'); hca.YLim = [0 7.999];
%hca = irf_panel('edi flux par'); hca.YLim = [0 7.999];
%hca = irf_panel('fpi flux par'); hca.YLim = [0 7.999];
%hca = irf_panel('fpi flux apar'); hca.YLim = [0 7.999];
%hca = irf_panel('fpi flux apar 2'); hca.YLim = [0 7.999];
%hca = irf_panel('fpi flux apar 12 comp'); hca.YLim = [0 3.999];
%hca = irf_panel('fpi edi diff'); hca.YLim = [-4 4];
%hca = irf_panel('fpi edi diff larger range'); hca.YLim = [-4 4];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];
%hca = irf_panel('E par dt'); hca.YLim = [-70 60];

doDoubleAxis = 0; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
  set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
  set(ax2,'YAxisLocation','right');
  set(ax2,'Color','none','box','off'); % color of axis      
  %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
  %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
  irf_timeaxis(ax2,'nolabels')
  ax2.XLabel.String = [];
  ax2.YLabel.String = {'\delta n/n'};
  ax2.YLabel.Interpreter = 'tex';    
  ax2.YTick = hca.YTick*1e-3/n0;  
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align(h)
%ax.YLabel.Position(1) = 1.07;

%set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required

%%
vmax = 90000e3;
nv = 5000;
v_vec = linspace(-vmax,vmax,nv);
dv = v_vec(2) - v_vec(1);

iff = 20;
[f0,params] = mms_20170706_135303.get_f0(iff);
n = params.n;
ntot = sum(n);
R = n(1)/ntot;
T = params.T;
vd = params.vd;
vt = params.vt;
n0 = sum(n);
ff4 = f0(v_vec,n,vd,vt);

colors = pic_colors('matlab');

hca = subplot(1,1,1);
plot(hca,v_vec*1e-6,ff4,'color',colors(3,:),'linewidth',1)
hca.XLim = [-25 25];
hca.XLabel.String = 'v_{||} (10^3 km/s)';
hca.YLabel.String = 'f_0 (s/m^4)';
hca.YLim = [0 2.6]*1e-3;

