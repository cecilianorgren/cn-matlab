% Data are obtained here:
% +table/separatrix_acceleration.m
% sep/compile_acc_pot_data.m

%% Get data from Cluster events, copied from +table/separatrix_acceleration.m
egedal.phi = [8 0.8 5 10 2 2 9 2.1 3 11 11 5 6 8 6 5 6 1]*1e3;
egedal.Te = [210 340 90 150 500 650 80 200 220 230 170 120 350 1500 70 400 150 100];
egedal.beta_lobe = [0.003 0.003 0.001 0.008 0.004 0.03 0.003 0.0003 0.001 0.0009 0.002 0.0003 0.0006 0.002 0.003 0.002 0.002 0.001];
egedal.beta_inflow = [0.008 0.21 0.026 0.018 0.22 0.15 0.011 0.038 0.034 0.054 0.030 0.0048 0.12 0.51 0.0094 0.064 0.011 0.017];

%% Get data from MMS events
localuser = datastore('local','user');
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

events = [1 3:20];
nevents = numel(events);
%acc_pot = nan(nevents,1);

ic = 1;
for ievent = 1:nevents
  event = events(ievent);
  sep.get_tints;
  load(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data')
  tsAccPot_nobg_abs = acc_pot_data.tsAccPot_nobg_abs;
  tsAccPot = tsAccPot_nobg_abs;
  
  for ic = 1:3
    time_tmp = tint_phi+0.1*[-1 1];
    acc_pot(ievent,ic) = max(tsAccPot{ic}.tlim(time_tmp).abs.data);
    time1{ievent,1} = time_tmp(1);
    time2{ievent,1} = time_tmp(2);
    beta_lobe(ievent,ic) = acc_pot_data.be_lobe;
    tepar_lobe(ievent,ic) = acc_pot_data.Tepar_lobe;
    tepar_sheet(ievent,ic) = acc_pot_data.Tepar_sheet;
    teperp_lobe(ievent,ic) = acc_pot_data.Teperp_lobe;
    teperp_sheet(ievent,ic) = acc_pot_data.Teperp_sheet;
    te_lobe(ievent,ic) = (acc_pot_data.Tepar_lobe + 2*acc_pot_data.Teperp_lobe)/3;
    te_sheet(ievent,ic) = (acc_pot_data.Tepar_sheet + 2*acc_pot_data.Teperp_sheet)/3;

    ne_lobe(ievent,ic) = acc_pot_data.n_lobe;
    ne_sep(ievent,ic) = acc_pot_data.n_sep;
    ne_sheet(ievent,ic) = acc_pot_data.n_sheet;
    
    B_lobe_(ievent,ic) = acc_pot_data.B_lobe;
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

clear event
event.event_id =     events;
event.time1 =        time1;
event.time2 =        time2;
event.phi =          acc_pot';
event.B =            B_lobe_';
event.n_lobe =       ne_lobe';
event.n_sheet =      ne_sheet';
event.n_sep =        ne_sep';
event.Tpar_lobe =    tepar_lobe';
event.Tpar_sheet =   tepar_sheet';
event.Tpar_sep =     [];
event.Tperp_lobe =   [];
event.Tperp_sheet =  [];
event.Tperp_sep =    [];
event.T_lobe =       te_lobe';
event.T_sheet =      te_sheet';
event.fce0 =         [];
event.fpe0 =         [];
event.beta_e_lobe =  beta_lobe';
event.beta_e_sheet = [];
event.beta_e_sep =   [];

% manual override some values:
event.T_sheet(nevents) = 1500;

if 0 % print out data as should be inserted in table in paper
  %%
  for ievent_ = 1:nevents
    ievent = ievent_;%events(ievent_);
    t1_utc = event.time1{ievent}.utc;
    t2_utc = event.time2{ievent}.utc;
    str = sprintf('%s - %s & %.0f & %.0f & %.0f & %.3f',t1_utc(1:23),t2_utc(15:23),100*round(event.phi(ievent)/100),10*round(event.T_lobe(ievent)/10),100*round(event.T_sheet(ievent)/100),event.beta_e_lobe(ievent));
    disp(str)
  end 
end

%% Figure for paper
figure(42)
nrows = 1;
ncols = 3;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

fontsize = 12;
isub = 1;
colors = mms_colors('matlab');

doEgedal = 1;
if 1 % phi vs beta e lobe
  hca = h(isub); isub = isub + 1;  
  hs1 = scatter(hca,event.beta_e_lobe,event.phi*1e-3,[]);
  hs1.CData = colors(1,:);
  if doEgedal
    hold(hca,'on')
    hs2 = scatter(hca,egedal.beta_lobe,egedal.phi*1e-3,[]);
    hold(hca,'off')  
    hs2.CData = colors(2,:);
    set(hca,'colororder',colors) 
    %irf_legend(hca,{'MMS';{'Cluster','(Egedal et al. 2015)'}},[0.98 0.98])
    %irf_legend(hca,{'MMS';{'Cluster'}},[0.98 0.98])
    
    legend(hca,{'MMS';'Cluster'},'Box','off')
  end
  %hl = legend(hca,{'MMS','Egedal 2015'});
  %hl.Box = 'off';   
  hca.XLabel.String =  '\beta_{e}^{lb}';
  hca.YLabel.String =  '\psi (keV)';
  %hca.YLabel.Interpreter = 'latex';
  hca.FontSize = fontsize;
  %hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
  hca.YLim(2) = ceil(max([event.phi egedal.phi]*1e-3)*1.01);
end
if 1 % phi/Telobe vs beta e lobe
  hca = h(isub); isub = isub + 1;
  hs1 = scatter(hca,event.beta_e_lobe,event.phi./(event.T_lobe),[]);
  
  hs1.CData = colors(1,:);
  if doEgedal
    hold(hca,'on')
    hs2 = scatter(hca,egedal.beta_lobe,egedal.phi./egedal.Te,[]);
    hold(hca,'off')  
    hs2.CData = colors(2,:);
    set(hca,'colororder',colors) 
    %irf_legend(hca,{'MMS';{'Cluster','(Egedal et al. 2015)'}},[0.98 0.98])
    %irf_legend(hca,{'MMS';{'Cluster'}},[0.98 0.98])
    legend(hca,{'MMS';'Cluster'},'Box','off')
  end
  %hca.XLabel.String =  '\beta_{e,lobe}';
  hca.XLabel.String =  '\beta_{e}^{lb}';
  hca.YLabel.String =  'e\psi/k_BT_{e}^{lb}';
  hca.FontSize = fontsize;
  %hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
  hca.YLim(2) = 120;
  hold(hca,'on')
  
  if 0
  phi_scaling = @(alpha,beta_power,nratio,beta) alpha.^2.*nratio./beta.^beta_power;
  beta_vec = linspace(0.0005,0.03,100);  
  alpha = [0.05 0.1 0.15 0.20]; % reconnection rate
  beta_power = 1;
  for ialpha = 1:numel(alpha)
    nratio = 3;  
    plot(hca,beta_vec,phi_scaling(alpha(ialpha),beta_power,nratio,beta_vec),'displayname',sprintf('R = %.2f, n_{in}/n_{out} = %.0f',alpha(ialpha),nratio),'color',colors(ialpha+2,:))
    nratio = 5;  
    plot(hca,beta_vec,phi_scaling(alpha(ialpha),beta_power,nratio,beta_vec),'displayname',sprintf('R = %.2f, n_{in}/n_{out} = %.0f',alpha(ialpha),nratio),'color',colors(ialpha+2,:),'linestyle','--')
  end
  end
  legend(hca)
  hold(hca,'off')
end
if 1 % phi/Tesheet vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,event.phi./(event.T_sheet),[])
  %hca.XLabel.String =  '\beta_{e,lobe}';
  hca.XLabel.String =  '\beta_{e}^{lb}';
  hca.YLabel.String =  'e\psi/k_BT_{e}^{sh}';
  set(hca,'colororder',colors) 
  %irf_legend(hca,{'MMS'},[0.98 0.98])  
  legend(hca,{'MMS'},'Box','off')
  hca.FontSize = fontsize;
  hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
  hca.XLim = h(2).XLim;
end

for ipanel = 1:npanels
  h(ipanel).XScale = 'lin';
  h(ipanel).XLim = [0 0.0301];
  h(ipanel).XTick = 0:0.01:0.05;
  h(ipanel).XGrid = 'on';
  h(ipanel).YGrid = 'on';
  %h(ipanel).XTick = 10.^[-3 -2 -1 0];
end
h(1).YTick = 0:2:12;
h(2).YTick = 0:20:120;

% event.phi = ;
% event.B = ;
% event.n_lobe = ;
% event.n_sheet = ;
% event.n_sep = ;
% event.Tpar_lobe = ;
% event.Tpar_sheet = ;
% event.Tpar_sep = ;
% event.Tperp_lobe = ;
% event.Tperp_sheet = ;
% event.Tperp_sep = ;
% event.fce0 = ;
% event.fpe0 = ;
% event.beta_e_lobe = ;
% event.beta_e_sheet = ;
% event.beta_e_sep = ;

%% 2nd figure for paper, discussion with scaling
figure(43)
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

fontsize = 12;
isub = 1;
colors = mms_colors('matlab');

doEgedal = 1;
if 1 % phi/Telobe vs beta e lobe
  hca = h(isub); isub = isub + 1;
  hs1 = scatter(hca,event.beta_e_lobe,event.phi./(event.T_lobe),[]);
  %hcb = colorbar('peer',hca); % not working for some reason
  
  hs1.CData = colors(1,:);
  if doEgedal
    hold(hca,'on')
    hs2 = scatter(hca,egedal.beta_lobe,egedal.phi./egedal.Te,[]);
    hold(hca,'off')  
    hs2.CData = colors(2,:);
    set(hca,'colororder',colors) 
    %irf_legend(hca,{'MMS';{'Cluster','(Egedal et al. 2015)'}},[0.98 0.98])
    %irf_legend(hca,{'MMS';{'Cluster'}},[0.98 0.98])
    legend(hca,{'MMS';'Cluster'},'Box','off');
  end
  %hca.XLabel.String =  '\beta_{e,lobe}';
  hca.XLabel.String =  '\beta_{e}^{lb}';
  hca.YLabel.String =  'e\psi/k_BT_{e}^{lb}';
  hca.FontSize = fontsize;
  %hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
  hca.YLim(2) = 120;
  hold(hca,'on')
  
  phi_scaling = @(alpha,beta_power,nratio,beta) alpha.^2.*nratio./beta.^beta_power;
  beta_vec = linspace(0.0005,0.03,100);  
  alpha = [0.05 0.1 0.15 0.20]; % reconnection rate
  beta_power = 1;
  for ialpha = 1:numel(alpha)
    nratio = 1/3;  
    plot(hca,beta_vec,phi_scaling(alpha(ialpha),beta_power,nratio,beta_vec),'displayname',sprintf('R = %.2f, n_{in}/n_{out} = %.2f',alpha(ialpha),nratio),'color',colors(ialpha+2,:))
    nratio = 1/5;  
    plot(hca,beta_vec,phi_scaling(alpha(ialpha),beta_power,nratio,beta_vec),'displayname',sprintf('R = %.2f, n_{in}/n_{out} = %.2f',alpha(ialpha),nratio),'color',colors(ialpha+2,:),'linestyle','--')
  end
  legend(hca);
  hold(hca,'off')
  hca.YTick = 0:20:120;
end
if 0 % nlobe vs nsheet
  hca = h(isub); isub = isub + 1;
  plot(hca,event.n_lobe,event.n_sheet,'o')
  hold(hca,'on')
  plot(hca,[0 hca.XLim(2)],[0 hca.YLim(2)],'k-')
  hold(hca,'off')
  hca.XLabel.String = 'n_{lb} (cm^{-3})';
  hca.YLabel.String = 'n_{sh} (cm^{-3})';
end
if 0
  hca = h(isub); isub = isub + 1;
  xval = event.beta_e_lobe;
  ycal = event.n_lobe./event.n_sheet;
  sval = event.phi./(event.T_lobe);
  cval = event.phi./(event.T_lobe);
  
  hs3 = scatter(hca,xval,ycal,sval*10,0*cval);  
  
  linestyles = {'-','--','-.',':'};
  phi_scaling = @(alpha,beta_power,nratio,beta) alpha.^2.*nratio./beta.^beta_power;
  beta_vec = linspace(0.0005,0.03,100);  
  nratio_vec = linspace(0.05,2,99);
  [Nratio,Beta] = meshgrid(nratio_vec,beta_vec);
  alpha = [0.1]; % reconnection rate
  beta_power = 1;
  hold(hca,'on')
  for ialpha = 1:numel(alpha)
    PHI = phi_scaling(alpha(ialpha),beta_power,Nratio,Beta);
    contour(hca,Beta,Nratio,PHI,[0:10:120],'linestyle',linestyles{ialpha}) 
  end
  hold(hca,'off')
end

for ipanel = 1:npanels
  %h(ipanel).XScale = 'lin';
  %h(ipanel).XLim = [0 0.0301];
  %h(ipanel).XTick = 0:0.01:0.05;
  h(ipanel).XGrid = 'on';
  h(ipanel).YGrid = 'on';
  h(ipanel).FontSize = 14;
  %h(ipanel).XTick = 10.^[-3 -2 -1 0];
end
%h(1).YTick = 0:2:12;

% event.phi = ;
% event.B = ;
% event.n_lobe = ;
% event.n_sheet = ;
% event.n_sep = ;
% event.Tpar_lobe = ;
% event.Tpar_sheet = ;
% event.Tpar_sep = ;
% event.Tperp_lobe = ;
% event.Tperp_sheet = ;
% event.Tperp_sep = ;
% event.fce0 = ;
% event.fpe0 = ;
% event.beta_e_lobe = ;
% event.beta_e_sheet = ;
% event.beta_e_sep = ;
