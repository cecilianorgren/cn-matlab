saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

events = 1:20;
nevents = numel(events);

%acc_pot = nan(nevents,1);
doPlot = 0;

ic = 1;
for ievent = 1:nevents
  event = events(ievent);
  sep.get_tints;
  load(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data')
  tsAccPot_nobg_abs = acc_pot_data.tsAccPot_nobg_abs;
  tsAccPot = tsAccPot_nobg_abs;
  tsAccPot_nobg_rel = acc_pot_data.tsAccPot_nobg_rel;
  acc_pot_all_nobg_abs{ievent} = tsAccPot_nobg_abs;
  acc_pot_all_nobg_rel{ievent} = tsAccPot_nobg_rel;
  
  for ic = 1:3
    [maxPot,iMaxPot] = max(tsAccPot{ic}.tlim(tint_phi+0.1*[-1 1]).abs.data);
    acc_pot(ievent,ic) = maxPot;
    acc_pot_abs(ievent,ic) = maxPot;
    
    [maxPotRel,iMaxPotRel] = max(tsAccPot_nobg_rel{ic}.tlim(tint_phi+0.1*[-1 1]).abs.data);
    acc_pot_rel(ievent,ic) = maxPotRel;
    
    iMaxPots(ievent,ic) = iMaxPot;
    beta_lobe(ievent,ic) = acc_pot_data.be_lobe;
    tepar_lobe(ievent,ic) = acc_pot_data.Tepar_lobe;
    tepar_sheet(ievent,ic) = acc_pot_data.Tepar_sheet;
    te_lobe(ievent,ic) = (acc_pot_data.Tepar_lobe + 2*acc_pot_data.Teperp_lobe)/3;
    te_sheet(ievent,ic) = (acc_pot_data.Tepar_lobe + 2*acc_pot_data.Teperp_lobe)/3;

    ne_lobe(ievent,ic) = acc_pot_data.n_lobe;
    ne_sep(ievent,ic) = acc_pot_data.n_sep;
    ne_sheet(ievent,ic) = acc_pot_data.n_sheet;
    
    B_lobe_(ievent,ic) = acc_pot_data.B_lobe;
  end
  if doPlot
    %%
    [h,h2] = initialize_combined_plot(2,1,1,0.6,'vertical');
    isub = 1;
    if 1 % timeseries of potential and point where they maximize, nobg
      hca = h(isub); isub = isub + 1;
      set(hca,'colororder',mms_colors('1234'))      
      celltmp = cellfun(@(x) x.tlim(tint_phi+0.1*[-1 1]).abs,acc_pot_all_nobg_abs{ievent},'UniformOutput',0);
      irf_plot(hca,celltmp,'comp')
      hold(hca,'on')
      %celltmppoint = cellfun(@(x) x(),celltmp,'UniformOutput',0);
      hlines = irf_plot(hca,{celltmp{1}(iMaxPot),celltmp{2}(iMaxPot),celltmp{3}(iMaxPot)},'comp');
      c_eval('hlines.Children(?).Marker = ''o'';',1:3)
      hold(hca,'off')
    end
    if 1 % timeseries of potential and point where they maximize
      hca = h(isub); isub = isub + 1;
      set(hca,'colororder',mms_colors('1234'))      
      celltmp = cellfun(@(x) x.tlim(tint_phi+0.1*[-1 1]).abs,acc_pot_all_nobg_rel{ievent},'UniformOutput',0);
      irf_plot(hca,celltmp,'comp')
      hold(hca,'on')
      %celltmppoint = cellfun(@(x) x(),celltmp,'UniformOutput',0);
      hlines = irf_plot(hca,{celltmp{1}(iMaxPot),celltmp{2}(iMaxPot),celltmp{3}(iMaxPot)},'comp');
      c_eval('hlines.Children(?).Marker = ''o'';',1:3)
      hold(hca,'off')
    end
    if 1 % plot shape of the 2D , not possible, need reduced for that.
    end
    drawnow
    pause(2)
  end
end

ic = 1;
acc_pot = acc_pot(:,ic);
beta_lobe = beta_lobe(:,ic);
B_lobe_ = B_lobe_(:,ic);
tepar_lobe = tepar_lobe(:,ic);
tepar_sheet = tepar_sheet(:,ic);
te_lobe = te_lobe(:,ic);
te_sheet = te_sheet(:,ic);
ne_lobe = ne_lobe(:,ic);
ne_sep = ne_sep(:,ic);
ne_sheet = ne_sheet(:,ic);

%% Plot for paper, including Cluster


%% Plot
nrows = 2;
ncols = 3;
npanels = ncols*nrows;
h = setup_subplots(nrows,ncols);
isub = 1;

%markersize = 240; markertype = '.';
markersize = 40;
markertype = 'o';
colormap(pic_colors('candy'))
%colormap('jet')
%colormap(pic_colors('c'))
color_data = events;

hca = h(isub); isub = isub + 1;
hscat = scatter(hca,beta_lobe,acc_pot*1e-3,markersize,color_data,'marker',markertype);
%hscat = scatter(hca,B_lobe_,acc_pot*1e-3,markersize,color_data,'marker',markertype);
hca.XLabel.String = '\beta_e^{lb}';
hca.YLabel.String = '\phi (keV)';

hca = h(isub); isub = isub + 1;
hscat = scatter(hca,beta_lobe,acc_pot./tepar_lobe,markersize,color_data,'marker',markertype);
hca.XLabel.String = '\beta_e^{lb}';
hca.YLabel.String = '\phi/T_e^{lb}';

hca = h(isub); isub = isub + 1;
hscat = scatter(hca,beta_lobe,acc_pot./tepar_sheet,markersize,color_data,'marker',markertype);
hca.XLabel.String = '\beta_e^{lb}';
hca.YLabel.String = '\phi/T_e^{sh}';

if 1 % ne_sheet vs ne_lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,ne_sep,ne_lobe,markersize,color_data,'marker',markertype)
  hca.XLabel.String = 'n_e^{sep}';
  hca.YLabel.String = 'n_e^{lb}';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0;
  hold(hca,'on')
  hpl = plot(hca,hca.XLim,hca.XLim,hca.XLim,hca.XLim*2,hca.XLim,hca.XLim*4);
  hold(hca,'off')
  %hleg = legend(hpl,{'n_e^{lb}/n_e^{sep} = 1','n_e^{lb}/n_e^{sep} = 2','n_e^{lb}/n_e^{sep} = 4'},'location','best');
  hleg = legend(hpl,{'1','2','4'},'location','best');
  hleg.Title.String = 'n_e^{lb}/n_e^{sep} = ';
end
if 1 % phi/Te_lobe vs ne_lobe/ne_sheet
  hca = h(isub); isub = isub + 1;
  scatter(hca,ne_lobe./ne_sep,acc_pot./tepar_lobe,markersize,color_data,'marker',markertype)
  hca.XLabel.String = 'n_e^{lobe}/n_e^{sep}';
  hca.YLabel.String = '\phi/T_e^{lb}';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0;
end
if 1 % ne_lobe/ne_sheet vs betae_lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,beta_lobe,ne_lobe./ne_sep,markersize,color_data,'marker',markertype)
  hca.XLabel.String = '\beta_e^{lb}';
  hca.YLabel.String = 'n_e^{lobe}/n_e^{sep}';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0;
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,beta_lobe,tepar_sheet./tepar_lobe,markersize,color_data,'marker',markertype)
  hca.XLabel.String = '\beta_e^{lb}';
  hca.YLabel.String = '\phi/T_e^{sh}';
end

for ipanel = 1:npanels
  h(ipanel).Box = 'on';
  h(ipanel).XGrid = 'on';
  h(ipanel).YGrid = 'on';
end


