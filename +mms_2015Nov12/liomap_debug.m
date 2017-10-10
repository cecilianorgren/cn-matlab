%% Compare single time distributions

time = times(1);

%% Plot flat skymap of all energy levels
tsPlot = tsFmap(1);
tsStop = tsZstop(1);
c_eval('gseBref = mean(gseB?.tlim(time+[-0.005 0.005]+-1).data,1);',ic)
c_eval('mvaBref = mean(mvaB?.tlim(time+[-0.005 0.005]+-1).data,1);',ic)
elevels = iE_;
energies = tsPlot.depend{1}(1,:);

fillval = 0;

tmpData = tsPlot.data;
tmpData(tmpData==fillval) = NaN;
tsPlot.data = tmpData;
npanels = numel(elevels);
nrows = 3;
ncols = 3;
npanels = ncols*nrows;

for isub = 1:npanels
  h(isub) = subplot(nrows,ncols,isub);
end
%% Ending position
for isub = 1:npanels
  hca = h(isub);
  pcolor(hca,squeeze(tsStop.depend{2}),squeeze(tsStop.depend{3}),squeeze(tsStop.data(1,elevels(isub),:,:))'*1e-3);
  hcb = colorbar('peer',hca);
end
%% Real data
for isub = 1:npanels
  hca = h(isub);
  mms.plot_skymap(hca,obsPDist.deflux,'tint',times(1),'flat','energy',energies(elevels(isub)),'vectors',{[-0.2984    4.3814    2.4994],'B'});
end
%% Lioville mapping data
for isub = 1:npanels
  hca = h(isub);
  mms.plot_skymap(hca,tsPlot.deflux,'tint',times(1),'flat','energy',energies(elevels(isub)),'vectors',{[-0.2984    4.3814    2.4994],'B'});
end

%% Comparison of real, liouville, and ending position
elevels = [12];

isub = 1;

% plottype = 'sphere';
% for ielevel = 1:numel(elevels)
%   hca = h(isub); isub = isub + 1;
%   mms.plot_skymap(hca,obsPDist.deflux,'tint',time(1)+-1,plottype,'energy',energies(elevels(ielevel)),'vectors',{[0 0.2 1],'B'});
% end
plottype = 'flat';
for ielevel = 1:numel(elevels)
  hca = h(isub); isub = isub + 1;
  mms.plot_skymap(hca,obsPDist.deflux,'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});
end
for ielevel = 1:numel(elevels)
  hca = h(isub); isub = isub + 1;
  mms.plot_skymap(hca,tsPlot.deflux,'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});  
end
for ielevel = 1:numel(elevels)
  hca = h(isub); isub = isub + 1;
  mms.plot_skymap(hca,tsZstop,'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});  
  hca.CLim = 1e3*[-30 30];
end
%%
for ielevel = 1:numel(elevels)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,squeeze(tsZstop.depend{2}),squeeze(tsZstop.depend{3}),squeeze(tsZstop.data(1,elevels(ielevel),:,:))'*1e-3);
  hcb = colorbar('peer',hca);
end

%% Figure plotting all the electron orbits.
tsPlot = tsFmap(1);
tsStop = tsZstop(1);
c_eval('gseBref = mean(gseB?.tlim(time+[-0.005 0.005]+-1).data,1);',ic)
c_eval('mvaBref = mean(mvaB?.tlim(time+[-0.005 0.005]+-1).data,1);',ic)
elevels = iE_;
energies = tsPlot.depend{1}(1,:);

% to get magnetic field xyz
mms_2015Nov12.Bmodel;

fillval = 0;

tmpData = tsPlot.data;
tmpData(tmpData==fillval) = NaN;
tsPlot.data = tmpData;
npanels = numel(elevels);
nrows = 3;
ncols = 3;
npanels = ncols*nrows;

for isub = 1:3;
  h(isub) = subplot(nrows,ncols,isub);
end
isub = isub + 1;
h(isub) = subplot(nrows,ncols,4:9);

% Comparison of real, liouville, and ending position
elevels = [12];
ielevel = 1;
isub = 1;
plottype = 'flat';

hca = h(isub); isub = isub + 1;
[ax,hcb1] =  mms.plot_skymap(hca,obsPDist.convertto('s^3/km^6'),'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});

hca = h(isub); isub = isub + 1;
[ax,hcb] =  mms.plot_skymap(hca,tsPlot.convertto('s^3/km^6'),'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});  
hca.CLim = h(isub-2).CLim;
hca.CLim = [0 800];

hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,tsZstop,'tint',times(1),plottype,'energy',energies(elevels(ielevel)),'vectors',{gseBref,'B'});  
hca.CLim = 1e3*[-30 30];

hca = h(isub); isub = isub + 1;
plot3(hca,xyz(:,1)*1e-3-40,(xyz(:,2)-xyz(end,2))*1e-3+30,xyz(:,3)*1e-3,'color',[0 0 0]);

h(1).CLim = h(2).CLim;

holdnow = 1;
for iEnergy = ielevel
  for iAzimuthal = 1:2:32
    for iPolar = 1:2:16
      if holdnow
        hold(hca,'on')
        holdnow = 0;
        hca.XLabel.String = 'L (km)';
        hca.YLabel.String = 'M (km)';
        hca.ZLabel.String = 'N (km)';
      end
      thisParticle = saveParticle{elevels(ielevel),iAzimuthal,iPolar};
      thisX = thisParticle.r(:,1)*1e-3; % L
      thisY = thisParticle.r(:,2)*1e-3; % M
      thisZ = thisParticle.r(:,3)*1e-3; % N
      if 1 % colormap according to f
        cmap = colormap(h(1)); % parula
        ncmap = size(cmap,1); % 64
        clim = h(1).CLim;
        colorind = ceil(ncmap*thisParticle.f*1e30/clim(end));
        if colorind == 0
          plotcolor = [1 1 1];
        elseif colorind > ncmap
          colorind = 64;
          plotcolor = cmap(fix(colorind),:);
        else
          plotcolor = cmap(fix(colorind),:);
        end
      else % colormap according to ending position
        if thisZ(end)<0;
          plotcolor = [0.1000    0.4000    1.0000];
        else
          plotcolor = [0.9500    0.7000         0];
        end
      end
      plot3(hca,thisX,thisY,thisZ,'color',plotcolor)
      pause(0.1)
    end
  end
end
hold(hca,'off')