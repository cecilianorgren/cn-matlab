function TS = find_acc_pot(fred_orig,varargin)
% FIND_PHI_ACC finds energy as a function fo time where f maximizes

[~,args,nargs] = axescheck(varargin{:});
     

% Default values
doVint = 0;
doEint = 0;
doResample = 0;
doRemoveAfter = 0;
doRelativeMinPeakProminence = 0;
%doMinPeakProminence = 0;
MinPeakProminence = 0;
MinPeakWidth = 0;
MaxPeakWidth = 0;

have_options = nargs > 1;
while have_options
  switch(lower(args{1}))
    case 'vint' % velocity interval to consider
      l = 2;
      vint = args{2};
      doVint = 1;
    case 'eint' % energy interval to consider
      l = 2;
      eint = args{2};
      doEint = 1;
    case 'removeafter'
      doRemoveAfter = 1;
    case 'minpeakprominence' % energy interval to consider
      l = 2;
      MinPeakProminence = args{2};
     % doMinPeakProminence = 1;
    case 'minpeakwidth' % energy interval to consider
      l = 2;
      MinPeakWidth = args{2};
      %doMinPeakWidth = 1;  
    case 'maxpeakwidth' % energy interval to consider
      l = 2;
      MaxPeakWidth = args{2};
      %doMinPeakWidth = 1;  
    case 'relativeminpeakprominence' % energy interval to consider
      l = 2;
      RelativeMinPeakProminence = args{2};      
      doRelativeMinPeakProminence = 1;
    case {'downsampling','downsample','resample'} % energy interval to consider
      l = 2;
      samplingstep = args{2};      
      doResample = 1;      
    otherwise
      fprintf('Unknown input argument: %s.\n',args{1})
      l = 1;
  end
  args = args((l+1):end);
  if isempty(args), break, end
end

if doResample
  fred = fred_orig.resample(fred_orig.time(1:samplingstep:end));
else
  fred = fred_orig;
end

% Perpare data
nt = fred.length;
phi = nan(nt,1);
f = fred.data;

units = irf_units;
velocity = fred.depend{1};
energy = sign(velocity)*0.5*units.me.*(velocity*1e3).^2/units.e;   

%if isempty(MinPeakWidth); MinPeakWidth = energy(end)-energy(1); end

vlim = [-Inf Inf];
if not(doRemoveAfter)
  if doEint
    f(abs(energy)<eint(1)) = 0;
    f(abs(energy)>eint(2)) = 0;
    vlim = sqrt(2*eint*units.e/units.me)*1e-3;
  end
  if doVint
    f(velocity<vint(1)) = 0;
    f(velocity>vint(2)) = 0;
    if vint(1)>vlim(1); vlim(1) = vint(1); end
    if vint(2)<vlim(2); vlim(2) = vint(2); end
  end
end

if 1 % findpeaks()
  for it = 1:nt
    if doRelativeMinPeakProminence
      MinPeakProminence = RelativeMinPeakProminence*max(f(it,:));
    end
    %if any([doRelativeMinPeakProminence doMinPeakProminence])
      [PKS,LOCS] = findpeaks(f(it,:),energy(it,:),'MinPeakProminence',MinPeakProminence,'MinPeakWidth',MinPeakWidth,'MaxPeakWidth',MaxPeakWidth);
    %else
    %  [PKS,LOCS] = findpeaks(f(it,:),energy(it,:));
    %end
    if not(isempty(LOCS))
      % find peak corresponding to highest energy           
      iloc = find(abs(LOCS)==max(abs(LOCS)));
      if numel(iloc)>1
        disp(sprintf('it = %g: Multiple peaks at same (plus/minus) energy. Using largest positive value: iloc(2). Consider setting vint.',it));
        iloc = iloc(2);
      end      
      phi(it) = LOCS(iloc);
      vel(it) = sign(LOCS(iloc))*sqrt(2*abs(LOCS(iloc))*units.e/units.me)*1e-3;      
    else
      phi(it) = NaN;
      vel(it) = NaN;
    end
  end
else % max()
  [val,ind] = max(f');
  phi = val;

  for it = 1:nt
    phi(it) = energy(ind(it));
    vel(it) = velocity(ind(it));
  end
end

if doRemoveAfter
  if doEint
    phi(phi<eint(1)) = NaN;
    phi(phi>eint(2)) = NaN;    
    vlim = sqrt(2*eint*units.e/units.me)*1e-3;
  end
  if doVint
    phi(vel<eint(1)) = NaN;
    phi(vel>eint(2)) = NaN;    
    if vint(1)>vlim(1); vlim(1) = vint(1); end
    if vint(2)<vlim(2); vlim(2) = vint(2); end
  end
end
vel(isnan(phi)) = NaN;

TS = irf.ts_scalar(fred.time,phi);
tsVel = irf.ts_scalar(fred.time,vel);

doPlot = 1;
if doPlot  
  %%  
  str_info = {};
  if doRelativeMinPeakProminence
    str_info{end+1,1} = sprintf('RelativeMinPeakProminence = %g',RelativeMinPeakProminence);
  else %if doMinPeakProminence
    str_info{end+1,1} = sprintf('MinPeakProminence = %g',MinPeakProminence);
  end
  if doResample
    str_info{end+1,1} = sprintf('Resampling to = (1:%g:end)',samplingstep);
  end
  
  
  h = irf_plot(3);  
  plotAccVel = 1;
  steps = [1 2 4 6];  % for downsampling
  
  if 1 % original fred
    hca = irf_panel('fred orig');  
    fred_min = 1e-6;
    fred_to_plot = fred_orig; 
    fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
    irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));    

    irf_timeaxis(hca);
    if plotAccVel % plot acc. vel
      hold(hca,'on')
      for istep = 1%:numel(steps)
        h_(istep) = irf_plot(hca,tsVel.resample(tsVel.time(1:steps(istep):end))*1e-3,'.-');
        legends{istep} = sprintf('1:%g:end',steps(istep));
      end    
      hold(hca,'off')
    end
    if 1 % plot vlim
      hold(hca,'on')
      irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(1)*[1 1])*1e-3,'--k');
      irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(2)*[1 1])*1e-3,'--k');
      hold(hca,'off')
    end
    hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
    if 0%plotAccVel
      hlegs = legend(h_,legends,'location','best');
      hlegs.Title.String = 'downsampled to:';
    end
  end
  
  if 1 % fred, potentially downsampled
    hca = irf_panel('fred');  
    fred_min = 1e-6;
    fred_to_plot = fred; 
    fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
    irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));    

    irf_timeaxis(hca);
    if plotAccVel % plot acc. vel
      hold(hca,'on')
      for istep = 1%:numel(steps)
        h_(istep) = irf_plot(hca,tsVel.resample(tsVel.time(1:steps(istep):end))*1e-3,'.-');
        legends{istep} = sprintf('1:%g:end',steps(istep));
      end    
      hold(hca,'off')
    end
    if 1 % plot vlim
      hold(hca,'on')
      irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(1)*[1 1])*1e-3,'--k');
      irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(2)*[1 1])*1e-3,'--k');
      hold(hca,'off')
    end
    hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
    if 0%plotAccVel
      hlegs = legend(h_,legends,'location','best');
      hlegs.Title.String = 'downsampled to:';
    end
    if 0 % exluded E/v interval is zero in distribution
      hca = irf_panel('fred 2');  
      fred_min = 1e-6;
      fred_to_plot = fred.clone(fred.time,f); 
      fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
      irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));    
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
      irf_timeaxis(hca);  
      if 1 % plot acc. vel
        hold(hca,'on')
        irf_plot(hca,tsVel*1e-3,'.-');
        hold(hca,'off')
      end
    end
    
    hleg = irf_legend(hca,str_info,[0.01 0.98],'k'); 
    c_eval('hleg(?).BackgroundColor = [1 1 1];',1:numel(hleg));
  end
  
  if 1 % plot acc. pot
    hca = irf_panel('acc. pot');  
    hold(hca,'on')
    for istep = 1:numel(steps)
      irf_plot(hca,TS.resample(TS.time(1:steps(istep):end)),'.-');
      legends{istep} = sprintf('1:%g:end',steps(istep));
    end
    hold(hca,'off')
    hca.YLabel.String = {'acc. pot.', '(eV)'};
    hlegs = legend(hca,legends,'location','best');
    hlegs.Title.String = {'additionally','downsampled to:'};
  end
  
  irf_plot_axis_align
  irf_zoom(h,'x',fred.time)
end

return