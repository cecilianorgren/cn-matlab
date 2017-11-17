ic = 1;
time_interval = irf.tint('2015-10-16T10:33:44.7/2015-10-16T10:33:45.8');
time_interval = irf.tint('2015-10-16T10:33:44.7/2015-10-16T10:33:44.8');
[time_indices,~] = dist.time.tlim(time_interval);

% Calculate pitch angle distribution
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)



%%
%mms.plot_int_projection(ePDist1,'t',t,'vlim',1e4,'xyz',[0,1,0;-1,0,0; 0,0,1],'vzint',2000*[-1,1],'nmc',100);
%mms.plot_projection(ePDist1,'t',t,'vlim',1e4,'xyz',[0,1,0;-1,0,0; 0,0,1],'vzint',2000*[-1,1],'nmc',100);

%% Plot single time particle distributions, 1 sc, 4 projections,
% Set up plot
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim = 15*1e3;
elevlim = 20;
strCMap = 'jet';
projclim = [0 5];  
 

for it = time_indices(10);time_indices(1):time_indices(end);
  time = dist.time(it);
  %if exist('hmark'); delete(hmark); end
  %hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  
  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  if 1 % Perpendicular plane    
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(1:23);
  end
  if 1 % B plane 1
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';
  end
  if 1 % B plane 2
    hca = h(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';
  end
  if 0 % Pitchangle distribution
    hca = h(isub); isub = isub + 1;
    plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
    hca.YScale = 'log'; hca.XScale = 'log';
    hca.YLabel.String = ['f_e (' ePitch.units ')'];
    hca.XLabel.String = 'E (eV)';
    hca.XLim = [10 4000];
    legend(hca,{'0','90','180'})
    hca.YTick = 10.^[-4:5];
    %hca.YLim = [1e0 1e5];
  end
  if 1 % Perpendicular plane, integrated distribution
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(1:23);
    hcb = colorbar('peer',hca);
  end
  if 1 % B plane 1, integrated distribution
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
    %mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';
    hcb = colorbar('peer',hca);
  end  
  if 1 % B plane 2
    hca = h(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
    %mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';
    hcb = colorbar('peer',hca);
  end
  if 0 % Pitchangle distribution
    hca = h(isub); isub = isub + 1;
    plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,:)));
    hca.YScale = 'log'; hca.XScale = 'log';
    hca.YLabel.String = ['f_e (' ePitch.units ')'];
    hca.XLabel.String = 'E (eV)';
    hca.XLim = [10 4000];
    labels = arrayfun(@(x) {num2str(x)},ePitch.depend{2});
    legend(hca,labels)
    hca.YTick = 10.^[-4:5];
    hca.YLim = [1e-4 1e5];
  end
  vzint = [8000 10000];
  if 1 % Perpendicular plane, integrated distribution
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vzint,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(1:23);
    ht = text(hca,hca.XLim(1),hca.YLim(2),sprintf('%g<v<%g (km/s)',vzint),'verticalalignment','top','horizontalalignment','left','BackgroundColor',[1 1 1]);
    hcb = colorbar('peer',hca);
  end
  if 1 % B plane 1, integrated distribution
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
    %mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vzint,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    ht = text(hca,hca.XLim(1),hca.YLim(2),sprintf('%g<v<%g (km/s)',vzint),'verticalalignment','top','horizontalalignment','left','BackgroundColor',[1 1 1]);
    hca.Title.String = '';
    hcb = colorbar('peer',hca);
  end  
  if 1 % B plane 2
    hca = h(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
    %mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vzint,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    ht = text(hca,hca.XLim(1),hca.YLim(2),sprintf('%g<v<%g (km/s)',vzint),'verticalalignment','top','horizontalalignment','left','BackgroundColor',[1 1 1]);
    hca.Title.String = '';
    hcb = colorbar('peer',hca);
  end
  for ii = 1:npanels
    colormap(h(ii),strCMap)
  end
  %cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',eventPath)
end

%% Plot single time particle distributions, 1 sc, many slices in vpar
% Set up plot
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim = 20*1e3;
elevlim = 20;
strCMap = 'jet';
projclim = [0 5];  
 

for it = time_indices(4);time_indices(1):time_indices(end);
  time = dist.time(it);
  %if exist('hmark'); delete(hmark); end
  %hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  
  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  Eint_vector = logspace(log10(3),3,5);
  vzint_vector = [sort(-sqrt(Eint_vector*units.eV*2/units.me)) sqrt(Eint_vector*units.eV*2/units.me)]*1e-3;
  %vzint_vector = -9000:2000:9000;
  
  for isubplot = 1:numel(vzint_vector)-1
  vzint = vzint_vector(isubplot+[0 1]);
  if 1 % Perpendicular plane, integrated distribution
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
    mms.plot_int_projection(hca,dist,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vzint,'scpot',scpot,'vlabel',vlabels,'colorbar',1);   
    hca.Title.String = '';
    hca.Title.String = timeUTC(1:23);
    ht = text(hca,hca.XLim(1),hca.YLim(2),sprintf('%.0f<v_{||}<%.0f (km/s)',vzint),'verticalalignment','top','horizontalalignment','left','BackgroundColor',[1 1 1]);
    %hcb = colorbar('peer',hca);
    %colormap(hca,strCMap)
  end
  end
  %cn.print(sprintf('ring_integrated_9slices_mms%g_%s',ic,time.utc),'path',eventPath)
end
