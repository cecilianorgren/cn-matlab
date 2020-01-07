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

results = cell(20,1);

%% Make/compile phi statistics
% All events have slightly different conditions for finding the most
% distinct beam. For example spacecraft, psd threshhold (10%/20%).
% All reduced distributions are saved to file, see get_acc_pot_new_rerun.m.
% Run through them all and just specify the conditions we want, and pick
% out the value.

events = 1:20;
nEvents = numel(events);
phi_all = zeros(nEvents,1);

for iEvent = 1:5%:nEvents
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
      cond.ic = 2;
      cond.tint_phi = tint_phi + 1.5*[-1 1];
      cond.dist = tlim(ef1D_nobg{cond.ic},cond.tint_phi);
      cond.flow = 0.2; % percentage of f_lobe  
      cond.minpeakheight = f_lobe*cond.flow;
      cond.minpeakprominence = 0.5*cond.minpeakheight;
      cond.vlow = 3000; 
      cond.eint = eint;      
    case 6
      continue
    case 7
      continue
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
      continue
    case 10
      continue
    case 11
      continue
    case 12
      continue
    case 13
      continue
    case 14
      continue
    case 15
      continue
    case 16
      continue
    case 17
      continue
    case 18
      continue
    case 19
      continue
    case 20
      continue
  end

  
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
  
  %% Find beam    
  [tmpTsAccPot,tmppPeakF] = find_acc_pot(cond.dist,...
    'eint',cond.eint,...
    'vint',cond.vint,...
    'minpeakheight',cond.minpeakheight,... % fraction of f_lobe
    'minpeakprominence',cond.minpeakprominence,... % smaller fraction of f_lobe
    'resample',cond.resample_timestep);
  % Collect data
  cond.ts_acc_pot = tmpTsAccPot;
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
  results{iEvent} = cond;
end
