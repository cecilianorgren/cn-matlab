% paper_fig_2d_vdfs
% Not completely selfcontained yet?
time = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
time = irf_time('2017-07-06T13:54:05.50Z','utc>EpochTT');
tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_psbl = irf.tint('2017-07-06T13:54:05.520Z/2017-07-06T13:54:05.640Z');
tint_cold = irf.tint('2017-07-06T13:54:13.000Z/2017-07-06T13:54:13.040Z')+0;

mms_id = 1;

%% Make reduced distribution, EH range
strTint = [irf_time(tint(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vint = [-Inf Inf];
vg = (-100:2:100)*1e3;
c_eval('eDist = ePDist?.tlim(tint_psbl);',mms_id)

% remove background
nSecondary = [5];
nPhoto = 0;
%[eDist_nobg] = mms.remove_edist_background(eDist_orig);
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''Nphotoe_art'',nPhoto,''ZeroNaN'',0);',1:numel(nSecondary))

c_eval('eDist_bgremoved = eDist_nobg?.tlim(tint_psbl);',mms_id)

vgi = [-800:50:800];
%c_eval('iDist = iPDist?.tlim(tint+[-0.05 +0.05]);',mms_id)

scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
zhat = irf.ts_vec_xyz(ePara.time,repmat([0 0 1],ePara.length,1));
ePerp1 = zhat.cross(ePara).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 40;
nMC = 500;
%tic; ef1D_ = eDist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
%tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
%tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
%tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

tic; ef1D_bgremoved = eDist_bgremoved.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2_bgremoved = eDist_bgremoved.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

% tic; if1D_ = iDist.reduce('1D',ePara,'vint',vint,'nMC',nMC); toc % reduced distribution along B
% tic; if2D_parperp1 = iDist.reduce('2D',ePara,ePerp1,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc 
% tic; if2D_parperp2 = iDist.reduce('2D',ePara,ePerp2,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc
% tic; if2D_perp1perp2 = iDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc

% Make pitch angle spectrograms
%ePitch = eDist.pitchangles(dmpaB1.resample(eDist),12);
%ePitch_bgremoved = eDist.pitchangles(dmpaB1.resample(eDist_bgremoved),12); 

%% Make reduced distribution, lobe/cold range
strTint_cold = [irf_time(tint_cold(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tint_cold(2),'epochtt>utc_HHMMSS')];
eint = [000 40000];
vint = [-Inf Inf];
vg = (-100:2:100)*1e3;
c_eval('eDist = ePDist?.tlim(tint_cold);',mms_id)
c_eval('eDist_bgremoved = eDist_nobg?.tlim(tint_cold);',mms_id)

vgi = [-800:50:800];
%c_eval('iDist = iPDist?.tlim(tint_cold+[-0.05 +0.05]);',mms_id)

scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
zhat = irf.ts_vec_xyz(ePara.time,repmat([0 0 1],ePara.length,1));
ePerp1 = zhat.cross(ePara).norm;
ePerp2 = ePara.cross(ePerp1).norm;

lowerelim = 40;
nMC = 500;
tic; ef1D_cold = eDist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_cold_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_cold_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_cold_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

tic; ef1D_cold_bgremoved = eDist_bgremoved.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_cold_parperp1_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
tic; ef2D_cold_parperp2_bgremoved = eDist_bgremoved.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
tic; ef2D_cold_perp1perp2_bgremoved = eDist_bgremoved.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

% tic; if1D_cold_ = iDist.reduce('1D',ePara,'vint',vint,'nMC',nMC); toc % reduced distribution along B
% tic; if2D_cold_parperp1 = iDist.reduce('2D',ePara,ePerp1,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc 
% tic; if2D_cold_parperp2 = iDist.reduce('2D',ePara,ePerp2,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc
% tic; if2D_cold_perp1perp2 = iDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'vg',vgi,'base','cart','nMC',nMC); toc


% Make pitch angle spectrograms
%ePitch_cold = eDist.pitchangles(dmpaB1.resample(eDist),12);
%ePitch_cold_bgremoved = eDist.pitchangles(dmpaB1.resample(eDist_bgremoved),12); 

%%
it = 2;
h = setup_subplots(1,3);
npanels = 3;
isub = 1;

if 1 % parperp1
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_parperp1_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
end
if 1 % par1perp2
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_parperp2_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
end
if 1 % perp1perp2
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_perp1perp2_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
end

links_2d = linkprop(h,{'XLim','YLim','CLim'});
for ipanel = 1:npanels
  h(ipanel).CLim = [-15 -9.5];
  axis(h(ipanel),'square');
  %h(ipanel).XLim = vg([1 end])*1e-3*1;
  %h(ipanel).YLim = vg([1 end])*1e-3*1;
  h(ipanel).XLim = 60*[-1 1];
  h(ipanel).YLim = 60*[-1 1];
end
colormap(pic_colors('candy4'))

for ipanel = 1:npanels
  h(ipanel).Position(1) = h(ipanel).Position(1)-0.05;
  axis(h(ipanel),'square');
  %h(ipanel).XLim = vg([1 end])*1e-3*1;
  %h(ipanel).YLim = vg([1 end])*1e-3*1;
  h(ipanel).XLim = 60*[-1 1];
  h(ipanel).YLim = 60*[-1 1];
end

hcb = findobj(gcf,'type','colorbar');
delete(hcb(3))
delete(hcb(2))
hcb(1).Position(1) = 0.9;
hcb(1).Position(2) = h(1).Position(2)+0.05;
hcb(1).Position(4) = h(1).Position(4)-0.05;
hcb(1).Position(3) = 0.02;

%% Plot, 2D + 1D, vphav overlaid
it = 2;
it = 2;
h = setup_subplots(1,3);

npanels = 3;
isub = 1;
vlim = [-70 70];
flim_2D = [-14.5 -9.5]; % log
flim_1D = [0 1.7]*1e-3;
flim_1D = [0 3.2]*1e-3;
manual = edi_event_manual_dt;
vph = [manual.vpar];
%vphmark = -9000; % vpar
vphmark = [min(vph),max(vph)]*1e-3;


% EDI parameters
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV
E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me)*1e-6; % Mm/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me)*1e-6; % Mm/s

colors = mms_colors('matlab');
edi_color = [0 0 0]; %colors(2,:);


if 1 % parperp1
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_parperp1_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)'; 
  h_all.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.CLim = flim_2D;
  position = hca.Position;
  delete(h_all.Colorbar)
  hca.Position = position;
  axis(hca,'square')
  hca.XTick = hca.YTick;
  if 1 % vph range
    hold(hca,'on')    
    for iv = 1:numel(vphmark)
      plot(hca,vphmark(iv)*[1 1],hca.YLim,'k')      
    end    
    patch(hca,[vphmark(1) vphmark(2) vphmark(2) vphmark(1) vphmark(1)],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              'b','facealpha',0.1)    
    hold(hca,'off')
  end
  if 1 % edi range
    hold(hca,'on')        
    %plot(hca,v_edi_minus*[1 1],hca.YLim,'k')
    %plot(hca,v_edi_plus*[1 1],hca.YLim,'k')          
    
    patch(hca,[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
    patch(hca,-[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
                
    hold(hca,'off')
  end
  set(hca,'children',flipud(get(hca,'children')))
  irf_legend(hca,{'towards';'X line'},[0.12 0.02],'color',[0 0 0],'fontsize',13,'horizontalalignment','center')
  irf_legend(hca,{'away from';'X line'},[0.85 0.02],'color',[0 0 0],'fontsize',13,'horizontalalignment','center')
end
if 1 % par1perp2
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_parperp2_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  h_all.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  h_all.Colorbar.YLabel.FontSize = 12;
  h_all.Colorbar.FontSize = 12;
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.CLim = flim_2D;
  axis(hca,'square')
  hca.XTick = hca.YTick;
  if 1
    hold(hca,'on')
    for iv = 1:numel(vphmark)
      plot(hca,vphmark(iv)*[1 1],hca.YLim,'k')
    end
    patch(hca,[vphmark(1) vphmark(2) vphmark(2) vphmark(1) vphmark(1)],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              'b','facealpha',0.1)    
    hold(hca,'off')
    hold(hca,'off')
  end
  if 1 % edi range
    hold(hca,'on')        
    %plot(hca,v_edi_minus*[1 1],hca.YLim,'k')
    %plot(hca,v_edi_plus*[1 1],hca.YLim,'k')          
    
    patch(hca,[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
    patch(hca,-[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
               
    hold(hca,'off')
  end
  set(hca,'children',flipud(get(hca,'children')))
end
if 1 % par
  hca = h(isub); isub = isub + 1;
  hhot = plot(hca,ef1D_bgremoved(it).depend{1}*1e-3,ef1D_bgremoved(it).data,'color',[0 0 0],'linewidth',1);
  hold(hca,'on')
  hcold = plot(hca,ef1D_cold_bgremoved(1).depend{1}*1e-3,ef1D_cold_bgremoved(1).data,'color',colors(1,:),'linewidth',1);
  hold(hca,'off')
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'f_e (s/m^4)';
  hca.XLim = vlim;  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = flim_1D;  
  %hleg = legend([hhot,hcold],{'time of EHs','lobe'});
  if 1
    hold(hca,'on')
    for iv = 1:numel(vphmark)
      plot(hca,vphmark(iv)*[1 1],hca.YLim,'k')
    end
    patch(hca,[vphmark(1) vphmark(2) vphmark(2) vphmark(1) vphmark(1)],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              'b','facealpha',0.1)    
    
    patch(hca,[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
    patch(hca,-[v_edi_minus v_edi_plus v_edi_plus v_edi_minus v_edi_minus],...
              [hca.YLim(1) hca.YLim(1) hca.YLim(2) hca.YLim(2) hca.YLim(1)],...
              edi_color,'facealpha',0.5)
            
    %set(hca,'children',flipud(get(hca,'children')))        
    hold(hca,'off')
  end
  %hleg = legend([hhot,hcold],{'time of EHs','lobe'});
end
if 0 % perp1perp2
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = ef2D_perp1perp2_bgremoved(it).plot_plane(hca);
  hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
end

h(3).XLim = [-35 35];
h(3).XLim = [-25 25];

legends = {'a)','b)','c)'};
legends_color = {'k','k','k'};
for ipanel = 1:3
  irf_legend(h(ipanel),legends{ipanel},[-0.15 0.99],'fontsize',15,'color',legends_color{ipanel});
end


h(1).Position(1) = 0.12;
h(2).Position(1) = 0.35;
h(3).Position(1) = 0.63;

h(1).Position(2) = 0.15;
h(2).Position(2) = 0.15;
h(3).Position(2) = 0.15;

h(1).Position(4) = 0.7;
h(2).Position(4) = 0.7;
h(3).Position(4) = 0.7;

for ipanel = 1:npanels
  h(ipanel).FontSize = 12;
end
colormap(pic_colors('candy4'))

%irf_legend(h(1),utc_warm,[0.02 0.75],'color',[0 0 0],'fontsize',13,'horizontalalignment','left')
%irf_legend(h(2),utc_warm,[0.98 0.98],'color',[0 0 0],'fontsize',13,'horizontalalignment','right')

% plot ring and annotate wamr/hot populations
if 1
  %%
v1warm = 15*sind(0:360)-6;
v2warm = 15*cosd(0:360)-0;
v1hot = 40*sind(0:360)+15;
v2hot = 40*cosd(0:360)+0;
hold(h(1),'on')
plot(h(1),v1warm,v2warm,'--k','linewidth',1)
plot(h(1),v1hot,v2hot,'--k','linewidth',1)
hold(h(1),'off')
hold(h(2),'on')
plot(h(2),v1warm,v2warm,'--k','linewidth',1)
plot(h(2),v1hot,v2hot,'--k','linewidth',1)
hold(h(2),'off')

annotation('textarrow',[0.405 0.425],[0.35 0.438],'string',{'warm/','heated','electrons'},'fontsize',13,'horizontalalignment','center');
annotation('textarrow',[0.405 0.42],[0.70 0.63],'string',{'hot','electrons'},'fontsize',13,'horizontalalignment','center');
end

%annotation('textarrow',[0.17 0.19],[0.78 0.78],'string',{'range of','obs. v_{ph}'},'fontsize',12);
annotation('textarrow',[0.175 0.19],[0.78 0.78],'string',{'antiparallel','EDI range'},'fontsize',13,'horizontalalignment','center');
annotation('textarrow',[0.246 0.225],[0.78 0.78],'string',{'parallel','EDI range'},'fontsize',13,'horizontalalignment','center');
annotation('textarrow',[0.195 0.195],[0.87 0.85],'string',{'range of observed v_{ph}'},'fontsize',13);

%annotation('textarrow',[0.175 0.185]+0.21,[0.50 0.50],'string',{'warm/heated','electrons'},'fontsize',13,'horizontalalignment','center');
%annotation('textarrow',[0.5 0.48],[0.35 0.45],'string',{'hot','electrons'},'fontsize',13,'horizontalalignment','center');

utc_warm = ef1D_(it).time.utc; utc_warm = utc_warm(12:22);
utc_cold = ef1D_cold(1).time.utc; utc_cold = utc_cold(12:22);
%annotation('textarrow',[0.68 0.69]+0.003,[0.45 0.39]+0.005,'string',{'warmer','population','observed','with EHs',['~' utc_warm]},'fontsize',13,'color',[0 0 0],'horizontalalignment','right','TextBackgroundColor',[1 1 1]);
%annotation('textarrow',[0.780 0.756]-0.005,[0.60 0.57]-0.2,'string',{['' utc_warm],'warm/heated','population','observed','with EHs'},'fontsize',13,'color',[0 0 0],'horizontalalignment','right','TextBackgroundColor',[1 1 1]);
annotation('textarrow',[0.777 0.765],[0.59 0.57]-0.2,'string',{['' utc_warm],'warm/heated','population','observed','with EHs'},'fontsize',13,'color',[0 0 0],'horizontalalignment','right','TextBackgroundColor',[1 1 1]);
%annotation('textarrow',[0.76 0.745],[0.60 0.57],'string',{'cold population','observed closer','to the lobe',['~' utc_cold]},'fontsize',13,'color',colors(1,:),'TextBackgroundColor',[1 1 1]);
annotation('textarrow',[0.764 0.748],[0.60 0.57]+0.1,'string',{['' utc_cold],'cold population','observed closer','to the lobe'},'fontsize',13,'color',colors(1,:),'TextBackgroundColor',[1 1 1]);
%%
links_2d = linkprop(h,{'XLim','YLim','CLim'});



for ipanel = 1:npanels
  h(ipanel).Position(1) = h(ipanel).Position(1)-0.05;
  axis(h(ipanel),'square');
  %h(ipanel).XLim = vg([1 end])*1e-3*1;
  %h(ipanel).YLim = vg([1 end])*1e-3*1;
  h(ipanel).XLim = 60*[-1 1];
  h(ipanel).YLim = 60*[-1 1];
end

hcb = findobj(gcf,'type','colorbar');
delete(hcb(3))
delete(hcb(2))
hcb(1).Position(1) = 0.9;
hcb(1).Position(2) = h(1).Position(2)+0.05;
hcb(1).Position(4) = h(1).Position(4)-0.05;
hcb(1).Position(3) = 0.02;