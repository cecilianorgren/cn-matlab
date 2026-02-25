function ax = plot_cross_section_psd(varargin)
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

[ax,args,nargs] = axescheck(varargin{:});

dist = args{1}; args = args(2:end); timeDist = dist.time;
B = args{1}; args = args(2:end); 

have_xlim = 0;
have_ylim = 0;
have_energies = 0;
plotScPot = 0;
tId = 1:dist.length;

if nargs > 1, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))    
    case 'tint'
      l = 2;
      tint = args{2};
      dist = dist.tlim(tint);       
      %[tId,~] = dist.time.tlim(tint); 
      if isempty(tId); irf.log('warning','No data for given time interval.'); return; end
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
    otherwise
      disp([''])
  end
  args = args(l+1:end);
  if isempty(args), break, end  
end

if plotScPot, 
  scPot = scPot.tlim(dist.time([1 end])); scPot = mean(scPot.data,1); 
  if isempty(scPot); 
    irf.log('warning','Empty spacecraft potential input.'); 
    plotScPot = 0;
  end  
end

B = B.resample(dist.time); 
if isempty(B); irf.log('warning','Empty magnetic field input.'); return; end
if isempty(dist); irf.log('warning','Empty distribution input.'); return; end

dist = squeeze(irf.nanmean(dist.data,1));
Bhat = mean(B.data,1)/mean(B.abs.data,1);

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
pospar = ones(length(phi),length(theta)); 
pospar(find(thetab > 15)) = NaN; % 2*15
posperp = ones(length(phi),length(theta)); 
posperp(find(thetab < 82.5)) = NaN; % 2*7.5
posperp(find(thetab > 97.5)) = NaN;
posapar = ones(length(phi),length(theta)); 
posapar(find(thetab < 165)) = NaN; % 2*15

distpar = dist;
distperp = dist;
distapar = dist;

for ii = 1:length(energy);
  distpar(ii,:,:)  = squeeze(dist(ii,:,:)).*pospar;
  distperp(ii,:,:) = squeeze(dist(ii,:,:)).*posperp;
  distapar(ii,:,:) = squeeze(dist(ii,:,:)).*posapar;
end

distpar =  squeeze(irf.nanmean(irf.nanmean(distpar,3),2))*1e30; %Convert units (s^3 cm^-6 -> s^3 km^-6)
distperp = squeeze(irf.nanmean(irf.nanmean(distperp,3),2))*1e30;
distapar = squeeze(irf.nanmean(irf.nanmean(distapar,3),2))*1e30;



% Plot psd
if isempty(ax), fig = figure; ax = axes; end

plot(ax,energy,distpar,'k',energy,distperp,'r',energy,distapar,'b');
ylabel(ax,'f_e (s^3 km^{-6})');
xlabel(ax,'E (eV)')
set(ax,'yscale','log');
set(ax,'xscale','log');
irf_legend(ax,{'0 deg'},[0.91 0.92],'color','k')
irf_legend(ax,{'90 deg'},[0.91 0.84],'color','r')
irf_legend(ax,{'180 deg'},[0.91 0.76],'color','b')

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
  plot(ax,scPot+[0 0],ax.YLim)
  hold(ax,'off')
  %if scPot<ax.XLim(1), 
    irf_legend(ax,['scPot =' num2str(scPot,'%.1f') ' V'],[0.05 0.05],'color','k'); 
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

  