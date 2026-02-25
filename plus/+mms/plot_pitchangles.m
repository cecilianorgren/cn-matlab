function ax = plot_pitchangles(varargin)
% MMS.PLOT_cross_section_psd Plots psd for pitchangles 0 90 180.
%
%  [ax,hcb] = MMS.PLOT_cross_section_psd(TSeriesDesDist,B,'Opt1',Arg1,...) - plot skymap.
%  MMS.PLOT_SKYMAP(AX,...) - plot in axes AX.
%
%  Options:
%    'tint' - plot data for time interval if tint.length = 2
%             or closest time if tint.length = 1
%    'energies' - [energy1,energy2,...] - mark given energies
%    'ylim' - [ymin ymax]
%    'xlim' - [xmin xmax]
%    'scpot' - spacecraft potential in TSeries format. Marks scPot f?r time 
%              interval of distribution.
%    'binsize' - pitch angle +- limit. for pitch angle = 90, and binsize =
%                7.5, particle in the range [82.5 97.5] are included (total)
%                angle is 15 deg

[ax,args,nargs] = axescheck(varargin{:});

dist = args{1}; args = args(2:end); timeDist = dist.time;
B = args{1}; args = args(2:end); 

pitchAngles = [0 90 180]; nPA = numel(pitchAngles);
pitchAnglesWidth = 7.5; % plus minus
have_xlim = 0;
have_ylim = 0;
have_energies = 0;
plotScPot = 0;
tId = 1:dist.length;
tIdB = 1:B.length;

if nargs > 1, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))    
    case 'tint'
      l = 2;
      tint = args{2};
      if tint.length == 1 % find closest time
        tId = find(abs(dist.time-tint) == min(abs(dist.time-tint)) );
        %tIdB = find(abs(B.time-tint) == min(abs(B.time-tint)) );
      else % find indices for given time interval
        [tId,~] = dist.time.tlim(tint); 
        if isempty(tId); irf.log('warning','No data for given time interval.'); return; end              
      end
      
      %if 0 % old
      %l = 2;
      %tint = args{2};
      %dist = dist.tlim(tint);       
      %[tId,~] = dist.time.tlim(tint); 
      %if isempty(tId); irf.log('warning','No data for given time interval.'); return; end
      %end
    case 'scpot'
      l = 2;
      scPot = args{2};  
      plotScPot = 1;
    case 'energies'
      l = 2;
      energies = args{2};  
      have_energies = 1;  
    case 'ylim'
      l = 2;
      ylim = args{2};  
      have_ylim = 1;  
    case 'xlim'
      l = 2;
      xlim = args{2};  
      have_xlim = 1; 
    case 'binsize'
      l = 2;
      pitchAnglesWidth = args{2};      
    case 'pitchangle'
      l = 2;
      pitchAngles = args{2};
      nPA = numel(pitchAngles);      
    otherwise
      disp([''])
  end
  args = args(l+1:end);
  if isempty(args), break, end  
end

if plotScPot, 
  if isa(scPot,'TSeries')
    if numel(tId) == 1
      % get the data from surrounding 30 ms
      scPot = scPot.tlim(dist.time(tId)+0.015*[-1 1]); 
      scPot = mean(scPot.data,1); 
    else
      scPot = scPot.tlim(dist.time(tId(1):tId(end))); 
      scPot = mean(scPot.data,1); 
    end
  elseif isnumeric(scPot) % assume direct value is given
    scPot = mean(scPot,1);
  end
  if isempty(scPot); 
    irf.log('warning','Empty spacecraft potential input.'); 
    plotScPot = 0;
  end  
end

if isa(B,'TSeries')
  if numel(tId) == 1 % one time given
    % get the data from surrounding 30 ms
    B = B.tlim(dist.time(tId)+0.015*[-1 1]);     
  else % take average of time interval given
    B = B.tlim(dist.time(tId(1):tId(end)));      
  end
  %irf_plot(B); return
  B = mean(B.data,1); 
elseif isnumeric(B) % assume direct value is given
  B = mean(B,1);
else
  irf.log('warning','Can''t recognize magnetic field format, it should be TSeries or numerical vector.'); return;
end
if isempty(B); irf.log('warning','Empty magnetic field input.'); return; end

%B = B.resample(dist.time); 

if isempty(dist); irf.log('warning','Empty distribution input.'); return; end

dist = squeeze(irf.nanmean(dist(tId).data,1));
Bhat = B/norm(B);

dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);

% Calculate angles between B and theta and phi
thetab = acosd(x*Bhat(1)+y*Bhat(2)+z*Bhat(3));

% Sort into bins
psdPA = cell(nPA,1);
for ii = 1:nPA
  % Put all 'positions' which are not within the right angle range to NaN
  positions = ones(length(phi),length(theta)); 
  positions(find(thetab > pitchAngles(ii) + pitchAnglesWidth)) = NaN;
  positions(find(thetab < pitchAngles(ii) - pitchAnglesWidth)) = NaN;
  distribution = dist;
  for jj = 1:length(energy);    
    % Multiply dist*positions, so all that are not within the right angle range becomes NaN
    distribution(jj,:,:)  = squeeze(dist(jj,:,:)).*positions;
  end
  distribution =  squeeze(irf.nanmean(irf.nanmean(distribution,3),2))*1e30; %Convert units (s^3 cm^-6 -> s^3 km^-6)
  psdPA{ii} = distribution;
end

colorOrder = [0 0 0; 1 0 0; 0 0 1; 0 0.8 0; 1 0.8 0.0; 1 0 1];
legPosX = 0.92;
% Plot psd
if isempty(ax), fig = figure; ax = axes; end

for ii = 1:nPA
  plot(ax,energy,psdPA{ii},'color',colorOrder(ii,:))
  if ii == 1; hold(ax,'on'); end
  if ii == nPA; hold(ax,'off'); end  
  irf_legend(ax,{[num2str(pitchAngles(ii)) ' deg']},[0.91 0.92-(ii-1)*0.08],'color',colorOrder(ii,:))
end

titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
title(ax,titleString);

ylabel(ax,'f_e (s^3 km^{-6})');
xlabel(ax,'E (eV)')
set(ax,'yscale','log');
set(ax,'xscale','log');
box(ax,'on')

ax.XTick = 10.^[-1:7];

if have_xlim;
  ax.XLim = xlim;
else
  ax.XLim = [10 5e3];
end

if have_ylim;
  ax.YLim = ylim;  
end

if plotScPot 
  hold(ax,'on')
  plot(ax,abs(scPot)+[0 0],ax.YLim,'color',[0.1 0.7 0.8])
  hold(ax,'off')
  %if scPot<ax.XLim(1), 
    irf_legend(ax,['scPot =' num2str(scPot,'%.1f') ' V'],[0.35 0.05],'color',[0.1 0.7 0.8]); 
  %end
end

% Plot energies
hold(ax,'on');
while have_energies
  plot(ax,energies(1)+[0 0],ax.YLim,'k--')
  energies = energies(2:end);  
  if isempty(energies), break, end  
end
hold(ax,'off');

  