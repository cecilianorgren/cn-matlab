%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

units = irf_units;
doSave = 1;
doPrint = 1;

events = 1:20;
nevents = numel(events);
loadAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential_new/'];
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential_new_rerun/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

allEvents = 1:nevents;
okEvents = [1 5 8];
badEvents = setdiff(allEvents,okEvents);

results = struct();

%% Make/compile phi statistics
% All events have slightly different conditions for finding the most
% distinct beam. For example spacecraft, psd threshhold (10%/20%).
% All reduced distributions are saved to file, see get_acc_pot_new_rerun.m.
% Run through them all and just specify the conditions we want, and pick
% out the value.

events = 1:20;
nEvents = numel(events);
phi_all = zeros(nEvents,1);

countEvent = 0;

for iEvent = [1 3:19]%:nEvents
  countEvent = countEvent + 1;
  event = events(iEvent); % set event id
  sep.get_tints; % get tints for that event id  
  % adjust tint_phi below, to get proepr interval
  
  %% Load saved acc data
  load(sprintf('%s/acc_pot_data_event_%g',loadAccPotPath,event),'acc_pot_data'); % acc_pot_data
  tint_dist = acc_pot_data.ef1D_orig{1}.time([1 end]);
  fields = fieldnames(acc_pot_data);
  for ifield = 1:numel(fields)
    eval(sprintf('%s = acc_pot_data.(''%s'');',fields{ifield},fields{ifield}))
  end
  
  Te_sheet = (Tepar_sheet + 2*Teperp_sheet)/3;
  vte_sheet = sqrt(2*units.eV*Te_sheet/units.me);
  
  Te_lobe = (Tepar_lobe + 2*Teperp_lobe)/3;  
  vte_lobe = sqrt(2*units.eV*Te_lobe/units.me);
  f_lobe = 1e6*n_lobe/sqrt(pi)/vte_lobe; % sm^-4
  
  PB = (B_lobe*1e-9).^2/2/units.mu0;
  PP = (n_lobe*1e6)*(units.eV*Te_lobe);
  beta_lobe = PP/PB;
  
  % Pick out reconnection quadrant: tail, Earth to the left
  quadrant = ones(2,2);  
  if vix > 0; quadrant(:,2) = 0; else quadrant(:,1) = 0; end      
  if bix > 0; quadrant(2,:) = 0; else quadrant(1,:) = 0; end
   
  eint = [000 40000];
  
  %% Select what parameters to go with
  switch iEvent
    case 1
      % all are ok
      cond.ic = 4;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.dist = ef1D_nobg{cond.ic};
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 3000; 
      cond.eint = eint;
    case 2
      % really small drift, consider removing
      cond.ic = 4;
      cond.dist = ef1D_nobg{cond.ic};
      cond.flow = 0.2; % 20% of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 3000; 
      cond.eint = eint;
    case 3
      % beam is nice (example event), but threshold of 20% is too high 
      cond.ic = 1;
      cond.tint_phi = tint_phi + [0 0.5];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.11; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 10000; 
      cond.eint = eint;
    case 4
      % beam is nice
      cond.ic = 1;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.18; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 10000; 
      cond.eint = eint;
    case 5
      % all are ok
      % 2 and 4 are slightly higher, but 1 and 3 are less wobbly
      cond.ic = 1;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 3000; 
      cond.eint = eint;      
    case 6
      cond.ic = 3;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 3000; 
      cond.eint = eint;   
    case 7
      cond.ic = 3;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.05; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 7000; 
      cond.eint = eint;   
    case 8
      cond.ic = 1;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;       
    case 9
      cond.ic = 4;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.11; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 10
      cond.ic = 1;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.18; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 11
      cond.ic = 1;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.25; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 12
      cond.ic = 1;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.10; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 13
      cond.ic = 2;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.10; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 14
      cond.ic = 2;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.10; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 15
      cond.ic = 2;
      cond.tint_phi = tint_phi + 0.0*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.15; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 16
      cond.ic = 1;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.15; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 17
      cond.ic = 3;
      cond.tint_phi = tint_phi + 0.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 18
      cond.ic = 4;
      cond.tint_phi = tint_phi + 0.1*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.1; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 19
      cond.ic = 3;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.1; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
    case 20 % consider removing due to low fbeam
      cond.ic = 4;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.1; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 5000; 
      cond.eint = eint;  
  end

  %cond.flow = 0.05; % percentage of f_lobe  
  %cond.minpeakheight = f_lobe*cond.flow;
  
  % set parallel velocity to incude in interval
  if any(quadrant([1 4])) % top left or bottom right, inflow is antiparallel to B
    vint = [-100000 -cond.vlow];
  else  % top right or bottom left, inflow is parallel to B
    vint = [cond.vlow 100000];
  end
  cond.vint = vint;
  cond.resample_timestep = 3;
  cond.vte_sheet = vte_sheet*1e-3;
  cond.vte_lobe = vte_lobe*1e-3;
  cond.f_lobe = f_lobe;
  cond.Te_lobe = Te_lobe;
  cond.Te_sheet = Te_sheet;
  cond.be_lobe = be_lobe;
  
  %% Find beam
  figure(21)
  [tmpTsAccPot,tmppPeakF,tmpTsFe] = find_acc_pot(cond.dist,...
    'eint',cond.eint,...
    'vint',cond.vint,...
    'minpeakheight',cond.minpeakheight,... % fraction of f_lobe
    'minpeakprominence',cond.minpeakprominence,... % smaller fraction of f_lobe
    'resample',cond.resample_timestep);
  % Collect data
  cond.ts_acc_pot = tmpTsAccPot;
  cond.ts_f_all = tmpTsFe;
  cond.peakF = tmppPeakF;
  % Adjust figure and print
  if 1
    climfraction = [0.01 1];
    clim = climfraction*f_lobe;
    hcf = gcf;
    hcf.Position = [0 100 700 800];
    h = findobj(hcf.Children,'type','axes');
    c_eval('irf_pl_mark(h(?),tint_phi,''r'');',1:numel(h))
    c_eval('h(?).CLim = log10(clim);',2:3)
    %cn.print(sprintf('event%g_tsAccPot%g_orig',event,iic),'path',printAccPotPath)    
  end
  
  %% Collect results
  results(countEvent).iEvent = iEvent;
  fields = fieldnames(cond);
  for ifield = 1:numel(fields)
    results(countEvent).(fields{ifield}) = cond.(fields{ifield});
  end
  pause
  
end

%% Scatter plot of phi vs f (all)
figure(22)
h = setup_subplots(1,1);
isub = 1;

cmap = colormap('jet');
cmap = pic_colors('candy');
colors = [interp1(1:size(cmap,1),cmap(:,1),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,2),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,3),linspace(1,64,18))'];

if 1 % f/flobe vs v/vbeam, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    xx = abs(results(ir).ts_acc_pot.data)/results(ir).Te_lobe;
    yy = results(ir).ts_f_all.data/results(ir).f_lobe;
        
    xx(abs(xx)==min(abs(xx))) = NaN;
    plot(hca,xx,yy,'o','color',colors(ir,:))
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = '\psi/T_{e,lobe}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
  %hca.XLim = [-2 2];
  hca.XGrid ='on';
  hca.YGrid ='on';
  %hca.XTick = [-3:0.5:3];
  hca.FontSize = 14;
  %hca.XLim = [0 5];
end

%% Plot all beams and statistics
ebeam = [cellfun(@(x) ( x.ebeam ), {results(:).peakF}, 'UniformOutput', false)]; ebeam = [ebeam{:}];
vbeam = [cellfun(@(x) ( x.vbeam ), {results(:).peakF}, 'UniformOutput', false)]; vbeam = [vbeam{:}];
fbeam = [cellfun(@(x) ( x.fbeam ), {results(:).peakF}, 'UniformOutput', false)]; fbeam = [fbeam{:}];

betalobe = [results.be_lobe];
flobe = [results.f_lobe];
vtelobe = [results.vte_lobe];
vtesheet = [results.vte_sheet];

cmap = colormap('parula');
colors = [interp1(1:size(cmap,1),cmap(:,1),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,2),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,3),linspace(1,64,18))'];
  
figure(22)
h = setup_subplots(3,3);
isub = 1;


if 1 % ebeam vs beta_lobe
  hca = h(isub); isub = isub + 1;
  plot(hca,betalobe,abs(ebeam),'o')
  hca.XLabel.String = '\beta_{lobe}';
  hca.YLabel.String = 'E_{beam}';
end
if 1 % ebeam vs beta_lobe
  hca = h(isub); isub = isub + 1;
  scatter(hca,betalobe,abs(ebeam),1e3*fbeam./flobe)
  hca.XLabel.String = '\beta_{lobe}';
  hca.YLabel.String = 'E_{beam}';
end
if 0 % fbeam/flobe vs vbeam/vtesheet
  hca = h(isub); isub = isub + 1;
  plot(hca,(abs(vbeam)./vtesheet)',(fbeam./flobe)','o')
  hca.XLabel.String = '|v_{beam}|/v_{te,sheet}';
  hca.YLabel.String = 'f_{beam}/f_{lobe}';
end
if 1 % f/flobe vs v/vtesheet, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      xx = results(ir).peakF.v;
    else
      xx = results(ir).peakF.v;
    end
    yy = results(ir).peakF.f_orig/results(ir).f_lobe;
    yy = results(ir).peakF.f_orig/fbeam(ir);
    xlim = 6500;
    yy(abs(xx)<xlim) = NaN;
    plot(hca,xx/vbeam(ir),yy)
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,sheet}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
end
if 1 % f/flobe vs v/vbeam, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      xx = results(ir).peakF.v;
    else
      xx = results(ir).peakF.v;
    end
    yy = results(ir).peakF.f_orig/results(ir).f_lobe;
    %yy = results(ir).peakF.f_orig/fbeam(ir);
    xlim = 6500;
    yy(abs(xx)<xlim) = NaN;
    plot(hca,xx/vbeam(ir),yy,'color',colors(ir,:))
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{beam}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
  hca.XLim = [-2 2];
end
if 1 % f/flobe vs v/vsheet, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      xx = -results(ir).peakF.v;
    else
      xx = results(ir).peakF.v;
    end
    yy = results(ir).peakF.f_orig/results(ir).f_lobe;
    %yy = results(ir).peakF.f_orig/fbeam(ir);
    xlim = 6500;
    yy(abs(xx)<xlim) = NaN;
    plot(hca,xx/vtesheet(ir),yy,'color',colors(ir,:))
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,sheet}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
  hca.XLim = [-2 2];
end
if 1 % f/flobe vs v/vlobe, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      xx = -results(ir).peakF.v;
    else
      xx = results(ir).peakF.v;
    end
    yy = results(ir).peakF.f_orig/results(ir).f_lobe;
    %yy = results(ir).peakF.f_orig/fbeam(ir);
    xlim = 6500;
    yy(abs(xx)<xlim) = NaN;
    plot(hca,xx/vtelobe(ir),yy,'color',colors(ir,:))
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,lobe}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
  hca.XLim = [-8 8];
end
if 1 % f/flobe vs v/vtesheet
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      plot(hca,-results(ir).peakF.v/results(ir).vte_sheet,results(ir).peakF.f_orig/results(ir).f_lobe)
    else
      plot(hca,results(ir).peakF.v/results(ir).vte_sheet,results(ir).peakF.f_orig/results(ir).f_lobe)
    end
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,sheet}';
  hca.YLabel.String = 'f/f_{lobe}';
  hca.YLim(2) = 0.4;
end
if 0 % f/flobe vs v/vbeam
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)    
    plot(hca,results(ir).peakF.v/vbeam(ir),results(ir).peakF.f_orig/results(ir).f_lobe)        
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,sheet}';
  hca.YLabel.String = 'f/f_{lobe}';
end
if 1 % f/flobe vs v/vbeam
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)        
    plot(hca,results(ir).peakF.v/vbeam(ir),results(ir).peakF.f_orig/fbeam(ir))    
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{beam}';
  hca.YLabel.String = 'f/f_{beam}';
end
if 0 % f/flobe vs E/Tsheet
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      semilogy(hca,sqrt(-results(ir).peakF.v/results(ir).vte_sheet),results(ir).peakF.f_orig/results(ir).f_lobe)
    else
      semilogy(hca,sqrt(results(ir).peakF.v/results(ir).vte_sheet),results(ir).peakF.f_orig/results(ir).f_lobe)
    end
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{te,sheet}';
  hca.YLabel.String = 'f/f_{lobe}';
end

%% For paper, statistics
clear event
% Get data from Cluster events, copied from +table/separatrix_acceleration.m
egedal.phi = [8 0.8 5 10 2 2 9 2.1 3 11 11 5 6 8 6 5 6 1]*1e3;
egedal.Te = [210 340 90 150 500 650 80 200 220 230 170 120 350 1500 70 400 150 100];
egedal.beta_lobe = [0.003 0.003 0.001 0.008 0.004 0.03 0.003 0.0003 0.001 0.0009 0.002 0.0003 0.0006 0.002 0.003 0.002 0.002 0.001];
egedal.beta_inflow = [0.008 0.21 0.026 0.018 0.22 0.15 0.011 0.038 0.034 0.054 0.030 0.0048 0.12 0.51 0.0094 0.064 0.011 0.017];

ebeam = [cellfun(@(x) ( x.ebeam ), {results(:).peakF}, 'UniformOutput', false)]; ebeam = [ebeam{:}];
vbeam = [cellfun(@(x) ( x.vbeam ), {results(:).peakF}, 'UniformOutput', false)]; vbeam = [vbeam{:}];
fbeam = [cellfun(@(x) ( x.fbeam ), {results(:).peakF}, 'UniformOutput', false)]; fbeam = [fbeam{:}];

betalobe = [results.be_lobe];
flobe = [results.f_lobe];
vtelobe = [results.vte_lobe];
vtesheet = [results.vte_sheet];

event.beta_e_lobe = betalobe;
event.phi = abs(ebeam);
event.T_lobe = [results.Te_lobe];
event.T_sheet = [results.Te_sheet];
for iev = 1:18
  event.time1{iev} = results(iev).tint_phi(1);
  event.time2{iev} = results(iev).tint_phi(2);
end

if 0 % print out data as should be inserted in table in paper
  %%
  for ievent = 1:18
    all_epoch_unix(ievent) = event.time1{ievent}.epochUnix;
  end
  [~,sort_ind] = sort(all_epoch_unix);
  all_events = 1:18;
  all_events = all_events(sort_ind);
  for ievent_ = 1:18
    ievent = all_events(ievent_);
    t1_utc = event.time1{ievent}.utc;
    t2_utc = event.time2{ievent}.utc;
    str = sprintf('%s - %s & %4.0f & %3.0f & %4.0f & %5.3f \\\\ \\hline',t1_utc(1:23),t2_utc(15:23),100*round(event.phi(ievent)/100),10*round(event.T_lobe(ievent)/10),100*round(event.T_sheet(ievent)/100),event.beta_e_lobe(ievent));
    disp(str)
  end 
end

figure(22)
h = setup_subplots(1,3);
npanels = 3;
isub = 1;

if 0 % f/flobe vs v/vbeam, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  %plot(hca,)
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
    
    legend(hca,{'MMS';'Cluster'},'Box','on')
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
    legend(hca,{'MMS';'Cluster'},'Box','on')
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
  legend(hca,{'MMS'},'Box','on')
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


irf_legend(h(1),'(a)',[-0.25 1.01],'fontsize',14);
irf_legend(h(2),'(b)',[-0.25 1.01],'fontsize',14);
irf_legend(h(3),'(c)',[-0.25 1.01],'fontsize',14);

%% For paper
ebeam = [cellfun(@(x) ( x.ebeam ), {results(:).peakF}, 'UniformOutput', false)]; ebeam = [ebeam{:}];
vbeam = [cellfun(@(x) ( x.vbeam ), {results(:).peakF}, 'UniformOutput', false)]; vbeam = [vbeam{:}];
fbeam = [cellfun(@(x) ( x.fbeam ), {results(:).peakF}, 'UniformOutput', false)]; fbeam = [fbeam{:}];

betalobe = [results.be_lobe];
flobe = [results.f_lobe];
vtelobe = [results.vte_lobe];
vtesheet = [results.vte_sheet];

cmap = colormap('parula');
colors = [interp1(1:size(cmap,1),cmap(:,1),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,2),linspace(1,64,18))',interp1(1:size(cmap,1),cmap(:,3),linspace(1,64,18))'];

figure(22)
h = setup_subplots(1,1);
isub = 1;

if 1 % f/flobe vs v/vbeam, remove all below 100 eV
  hca = h(isub); isub = isub + 1;
  for ir = 1:numel(results)
    if results(ir).peakF.vbeam<0
      xx = results(ir).peakF.v;
    else
      xx = results(ir).peakF.v;
    end
    yy = results(ir).peakF.f_orig/results(ir).f_lobe;
    %yy = results(ir).peakF.f_orig/fbeam(ir);
    xlim = 6500;
    yy(abs(xx)<xlim) = NaN;
    plot(hca,xx/vbeam(ir),yy,'color',colors(ir,:),'linewidth',1.5)
    %plot(hca,xx/vbeam(ir),yy)
    if ir == 1
      hold(hca,'on')
    end
  end
  hold(hca,'off')
  hca.XLabel.String = 'v/v_{beam}';
  hca.YLabel.String = 'f/f_{lobe}';
  %hca.YLim(2) = 0.4;
  hca.XLim = [-2 2];
  hca.XGrid ='on';
  hca.YGrid ='on';
  hca.XTick = [-3:0.5:3];
  hca.FontSize = 14;
end
