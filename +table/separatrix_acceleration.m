events = cell(1,1);


event.phi =          [4105,  1376,  7394,  2072,  10234, 5504,  3193,  2439,  1637,  1777,  1137,  1922,  1087,  1637,  1504,  921,   1376,  1777];
event.B =            [16,    23,    16,    16,    16,    16,    19,    19,    19,    18,    18,    18,    19,    23,    23,    23,    23,    23];
event.n_lobe =       [0.015, 0.045, 0.015, 0.021, 0.021, 0.021, 0.044, 0.055, 0.055, 0.079, 0.079, 0.079, 0.031, 0.059, 0.059, 0.045, 0.045, 0.045];
event.n_sheet =      [0.031, 0.033, 0.016, 0.028, 0.020, 0.028, 0.112, 0.112, 0.057, 0.050, 0.050, 0.050, 0.171, 0.112, 0.112, 0.033, 0.033, 0.033];
event.n_sep =        [0.004, 0.018, 0.009, 0.010, 0.006, 0.010, 0.013, 0.023, 0.022, 0.062, 0.049, 0.017, 0.020, 0.012, 0.012, 0.012, 0.019, 0.012];
event.Tpar_lobe =    [295,   160,   295,   189,   189,   189,   231,   397,   397,   235,   235,   235,   396,   105,   105,   160,   160,   160];
event.Tpar_sheet =   [3337,  1180,  2936,  4437,  3452,  4437,  4210,  4210,  3859,  3620,  3620,  3620,  1223,  1307,  1307,  1180,  1180,  1180];
event.Tpar_sep =     [2345,  1211,  3980,  3877,  5173,  4257,  2472,  2559,  2290,  2015,  2342,  1086,  1158,  1413,  1413,  1099,  951,   1153];
event.Tperp_lobe =   [253,   140,   253,   216,   216,   216,   239,   290,   290,   213,   213,   213,   327,   113,   113,   140,   140,   140];
event.Tperp_sheet =  [3538,  995,   2626,  4469,  3105,  4469,  3381,  3381,  3555,  3262,  3262,  3262,  1052,  1073,  1073,  995,   995,   995];
event.Tperp_sep =    [1889,  991,   2945,  3258,  3050,  3492,  1865,  2153,  1930,  1379,  1659,  726,   952,   944,   944,   784,   620,   709];
event.fce0 =         [458,   633,   458,   457,   457,   457,   542,   541,   541,   517,   517,   517,   532,   642,   642,   633,   633,   633];
event.fpe0 =         [1585,  1619,  1145,  1501,  1269,  1501,  2999,  2999,  2144,  2004,  2004,  2004,  3708,  3002,  3002,  1619,  1619,  1619];
event.beta_e_lobe =  [0.006, 0.005, 0.006, 0.007, 0.007, 0.007, 0.011, 0.019, 0.019, 0.021, 0.021, 0.021, 0.012, 0.005, 0.005, 0.005, 0.005, 0.005];
event.beta_e_sheet = [0.188, 0.027, 0.069, 0.208, 0.110, 0.208, 0.216, 0.216, 0.286, 0.276, 0.276, 0.276, 0.227, 0.094, 0.094, 0.027, 0.027, 0.027];
event.beta_e_sep =   [0.014, 0.016, 0.048, 0.048, 0.032, 0.061, 0.021, 0.046, 0.042, 0.104, 0.111, 0.013, 0.018, 0.010, 0.010, 0.008, 0.010, 0.008];

egedal.phi = [8 0.8 5 10 2 2 9 2.1 3 11 11 5 6 8 6 5 6 1]*1e3;
egedal.Te = [210 340 90 150 500 650 80 200 220 230 170 120 350 1500 70 400 150 100];
egedal.beta_lobe = [0.003 0.003 0.001 0.008 0.004 0.03 0.003 0.0003 0.001 0.0009 0.002 0.0003 0.0006 0.002 0.003 0.002 0.002 0.001];
egedal.beta_inflow = [0.008 0.21 0.026 0.018 0.22 0.15 0.011 0.038 0.034 0.054 0.030 0.0048 0.12 0.51 0.0094 0.064 0.011 0.017];
%% Plot
n_events = numel(event.phi);

figure(41)

nrows = 3;
ncols = 3;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

if 1 % phi vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,event.phi,[],1:n_events)
  hca.XLabel.String =  'beta e lobe';
  hca.YLabel.String =  'phi';
end
if 1 % phi vs beta inflow, incl egedal
  hca = h(isub); isub = isub + 1;
  scatter(hca,egedal.beta_inflow,egedal.phi,[],'r')    
  hca.XLabel.String =  'beta inflow';
  hca.YLabel.String =  'phi';
  hca.Title.String = 'Egedal 2015';
end
if 1 % phi vs beta e inflow, egedal
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,event.phi,[],'k')
  hold(hca,'on')
  scatter(hca,egedal.beta_lobe,egedal.phi,[],'r')
  hold(hca,'off')  
  %hl = legend(hca,{'MMS','Egedal 2015'});
  %hl.Box = 'off';  
  set(hca,'colororder',[0 0 0; 1 0 0])
  irf_legend(hca,{'MMS';'Egedal 2015'},[0.98 0.98])
  hca.XLabel.String =  'beta e lobe';
  hca.YLabel.String =  'phi';
  
end
if 0 % phi vs fce0/fpe0
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.fce0./event.fpe0,event.phi,[],1:n_events)
  hca.XLabel.String =  'fce0/fpe0';
  hca.YLabel.String =  'phi';
end
if 1 % phi/Telobe vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,3*event.phi./(event.Tpar_lobe + 2*event.Tperp_lobe),[],1:n_events)
  hca.XLabel.String =  'beta e lobe';
  hca.YLabel.String =  'phi/Te lobe';
end
if 0 % phi/Telobe vs fce0/fpe0
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.fce0./event.fpe0,event.phi,[],1:n_events)
  hca.XLabel.String =  'fce0/fpe0';
  hca.YLabel.String =  'phi/Te lobe';
end
if 1 % phi/Tesheet vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,3*event.phi./(event.Tpar_sheet + 2*event.Tperp_sheet),[],1:n_events)
  hca.XLabel.String =  'beta e lobe';
  hca.YLabel.String =  'phi/Te sheet';
end
if 0 % fce0/fpe0 vs beta e lobe 
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.fce0./event.fpe0,event.beta_e_lobe,[],1:n_events)
  hca.XLabel.String =  'beta e lobe';
  hca.YLabel.String =  'fce0/fpe0';
end
if 1 % phi/Telobe vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,3*(event.Tpar_lobe + 2*event.Tperp_lobe),event.phi,[],1:n_events)
  hca.XLabel.String =  'T e sheet';
  hca.YLabel.String =  'phi';  
end
%% Figure for presentation
figure(42)
nrows = 3;
ncols = 1;
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
  hs1 = scatter(hca,event.beta_e_lobe,event.phi,[]);
  hs1.CData = colors(1,:);
  if doEgedal
    hold(hca,'on')
    hs2 = scatter(hca,egedal.beta_lobe,egedal.phi,[]);
    hold(hca,'off')  
    hs2.CData = colors(2,:);
    set(hca,'colororder',colors) 
    irf_legend(hca,{'MMS';{'Cluster','(Egedal et al. 2015)'}},[0.98 0.98])
  end
  %hl = legend(hca,{'MMS','Egedal 2015'});
  %hl.Box = 'off';   
  hca.XLabel.String =  '\beta_{e,lobe}';
  hca.YLabel.String =  '\phi (keV)';
  hca.FontSize = fontsize;
  %hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
end
if 1 % phi/Telobe vs beta e lobe
  hca = h(isub); isub = isub + 1;
  hs1 = scatter(hca,event.beta_e_lobe,event.phi./((event.Tpar_lobe + 2*event.Tperp_lobe)/3),[]);
  
  hs1.CData = colors(1,:);
  if doEgedal
    hold(hca,'on')
    hs2 = scatter(hca,egedal.beta_lobe,egedal.phi./egedal.Te,[]);
    hold(hca,'off')  
    hs2.CData = colors(2,:);
    set(hca,'colororder',colors) 
    irf_legend(hca,{'MMS';{'Cluster','(Egedal et al. 2015)'}},[0.98 0.98])
  end
  hca.XLabel.String =  '\beta_{e,lobe}';
  hca.YLabel.String =  '\phi/T_{e,lobe}';
  hca.FontSize = fontsize;
  %hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
end
if 1 % phi/Tesheet vs beta e lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,event.beta_e_lobe,event.phi./((event.Tpar_sheet + 2*event.Tperp_sheet)/3),[])
  hca.XLabel.String =  '\beta_{e,lobe}';
  hca.YLabel.String =  '\phi/T_{e,sheet}';
  set(hca,'colororder',colors) 
  irf_legend(hca,{'MMS'},[0.98 0.98])
  hca.FontSize = fontsize;
  hca.XLim(1) = min([event.beta_e_lobe, 0.004]);
  hca.Box = 'on';
  hca.XLim = h(2).XLim;
end

for ipanel = 1:npanels
  h(ipanel).XScale = 'lin';
  %h(ipanel).XTick = 10.^[-3 -2 -1 0];
end

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
