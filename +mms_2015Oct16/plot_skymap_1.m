tintTimeplot = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintTimeplot = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:35.00Z');
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');

tint = irf.tint('2015-10-16T10:33:30.40Z',0.04)+(-0.17-0.1);
sc = 2;

% Downsample ion moments
c_eval('fs = 1/(ni?brst.time(2)-ni?brst.time(1));',sc); fny = fs/2;
c_eval('pi?_lowres = irf_filt(Pi?brst,0,fny/2,fs,5);',sc)
c_eval('Ti?_lowres = irf_filt(Ti?brst,0,fny/2,fs,5);',sc)
c_eval('ni?_lowres = irf_filt(ni?brst,0,fny/2,fs,5);',sc)

% Plot a skymap fxor one time and several energy levels with vph and B0 inserted

% Get vectors for a small time interval
% I didn't find a mean function for TSeries, like 
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.time).tlim(tint).data);',sc); hatB0 = double(irf_norm(B0));
c_eval('E0 = mean(dslE?brst.tlim(tint).data);',sc); hatE0 = double(irf_norm(E0));
c_eval('vi = mean(vi?brst.resample(dslE?brst.time).tlim(tint).data);',sc); hatVi = double(irf_norm(vi));
c_eval('ve = mean(ve?brst.resample(dslE?brst.time).tlim(tint).data);',sc); hatVe = double(irf_norm(ve));
c_eval('vExB = mean(vExB?brst.resample(dslE?brst.time).tlim(tint).data);',sc); hatVExB = double(irf_norm(vExB));

% Set up coordinates for skymap plot
c_eval('desDist = desDist?;',sc);
c_eval('disDist = disDist?;',sc);

% Set up coordinate system separately for ions and electrons if they for
% some reason would be different from each other.
r = 1; % radius of sphere
% Electrons
phi_edges = linspace(0,2*pi,size(desDist.data,3)+1); % azimuthal angle bin edges?
theta_edges = linspace(0,pi,size(desDist.data,4)+1); % polar angle bin edges?
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
eX = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles were coming from
eY = -r*sin(THETA).*sin(PHI);
eZ = -r*cos(THETA);
% Ions
phi_edges = linspace(0,2*pi,size(disDist.data,3)+1); % azimuthal angle bin edges?
theta_edges = linspace(0,pi,size(disDist.data,4)+1); % polar angle bin edges?
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
iX = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles were coming from
iY = -r*sin(THETA).*sin(PHI);
iZ = -r*cos(THETA);

% Plot a few energy levels on a skymap for both ion and electron data
nrows = 3;
ncols = 4;
nSubPlots = nrows*ncols;
for k = 1:nSubPlots; h(k) = subplot(nrows,ncols,k); end

% Choose which energy levels to plot
ionEnergyLevels = [10 15 20 25]; % between 1:32
electronEnergyLevels = [8 11 12 13]; % between 1:32

[etId,~] = desDist.time.tlim(tint); 
[itId,~] = disDist.time.tlim(tint+0.1*[-1 1]); 


% Plot ion and electron moments
hca = subplot(nrows,ncols,1:4);
axes(hca);
c_eval('veline = irf_plot(hca,ve?brst);',sc); hold(hca,'on'); 
c_eval('viline = irf_plot(hca,vi?brst);',sc); hold(hca,'off');
%irf_legend(hca,{'v_{ex}','v_{ey}','v_{ez}','v_{ix}','v_{iy}','v_{iz}'},[0.9 0.9])
legend([veline(:); viline(:)],{'v_{ex}','v_{ey}','v_{ez}','v_{ix}','v_{iy}','v_{iz}'},'location','EastOutside')
irf_pl_mark(hca,tint.epochUnix','red')
irf_zoom(hca,'x',[vi3.time.start.epochUnix vi3.time.stop.epochUnix])
title(hca,'Ion and electron velocity moments')
irf_zoom(hca,'x',tintTimeplot)
hca.YLabel.String = 'v [km/s]';

hcaB = axes;
hcaB.Position = hca.Position;
hcaB.Position(2) = hcaB.Position(2)+hcaB.Position(4)*0.5;
hcaB.Position(4) = hcaB.Position(4)*0.5;
hca.Position(4) = hca.Position(4)*0.5;

if 1 % B
  c_eval('irf_plot(hcaB,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?brst.abs},''comp'');',sc)
  hcaB.YLabel.String = 'B_{DMPA} [nT]';
  irf_legend(hcaB,{'B_x','B_y','B_z','|B|'},[0.98 0.95])  
else % E  
  c_eval('irf_plot(hcaB,{dslE?brst.x,dslE?brst.y},''comp'');',ic)
  hcaB.YLabel.String = {'E_{DSL}','[mV/m]'};
  irf_legend(hcaB,{'E_x','E_y'},[0.98 0.95]); 
end
irf_zoom(hcaB,'x',tintTimeplot)
hcaB.XLabel = [];
hcaB.XTickLabel = [];

% Skymap distributions
for k = 1:ncols    
    % Plot electron skymap data in second row
    hca = h(ncols*1+k);    
    axes(hca)
    C = squeeze(nanmean(desDist.data(etId,electronEnergyLevels(k),:,:),1))';
    hs = surf(hca,eX,eY,eZ,C);
    hc = colorbar('peer',hca);
    axis(hca,'square')
    axis(hca,'equal')
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    %titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(electronEnergyLevels(k))]};
    titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS') ' + ' num2str(tint.stop-tint.start) ' s'],['Energy level = ' num2str(electronEnergyLevels(k))]};
    hca.Title.String = titleString;   
    hc.YLabel.String = 'Electron phase spase density';
    shading flat;
    
    % Plot vectors
    hold(hca,'on');
    hold on
    scale = 1.5;    
    quiver3(hca,0,0,0,hatVExB(1),hatVExB(2),hatVExB(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVi(1),hatVi(2),hatVi(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVe(1),hatVe(2),hatVe(3),scale,'linewidth',2)
    quiver3(hca,-scale*hatB0(1),-scale*hatB0(2),-scale*hatB0(3),hatB0(1),hatB0(2),hatB0(3),2*scale,'linewidth',2)%quiver3(hca,0,0,0,hatB0(1),hatB0(2),hatB0(3),scale,'linewidth',2)
    quiver3(hca,-scale*hatE0(1),-scale*hatE0(2),-scale*hatE0(3),hatE0(1),hatE0(2),hatE0(3),2*scale,'linewidth',2)
    
    % Label vectors
    scale = 1.7;
    text(scale*hatVExB(1),scale*hatVExB(2),scale*hatVExB(3),'v_{ExB}','fontsize',14)
    text(scale*hatVi(1),scale*hatVi(2),scale*hatVi(3),'v_{i}','fontsize',14)
    text(scale*hatVe(1),scale*hatVe(2),scale*hatVe(3),'v_{e}','fontsize',14)
    text(scale*hatB0(1),scale*hatB0(2),scale*hatB0(3),'B_{0}','fontsize',14)
    text(scale*hatE0(1),scale*hatE0(2),scale*hatE0(3),'E_{0}','fontsize',14)
    hold(hca,'off');
    hold off
    
    % Plot ion skymap data in bottom row
    hca = h(ncols*2+k);    
    axes(hca)
    C = squeeze(nanmean(disDist.data(itId,ionEnergyLevels(k),:,:),1))';
    hs = surf(hca,iX,iY,iZ,C);
    hc = colorbar('peer',hca);
    axis(hca,'square')
    axis(hca,'equal')
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    %titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(ionEnergyLevels(k))]};
    titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS') ' + ' num2str(tint.stop-tint.start) ' s'],['Energy level = ' num2str(ionEnergyLevels(k))]};
    hca.Title.String = titleString;   
    hc.YLabel.String = 'Ion phase spase density';
    shading flat;
    
    % Plot vectors
    hold(hca,'on');
    hold on
    scale = 1.5;
    
    quiver3(hca,0,0,0,hatVExB(1),hatVExB(2),hatVExB(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVi(1),hatVi(2),hatVi(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVe(1),hatVe(2),hatVe(3),scale,'linewidth',2)        
    quiver3(hca,-scale*hatB0(1),-scale*hatB0(2),-scale*hatB0(3),hatB0(1),hatB0(2),hatB0(3),2*scale,'linewidth',2) %quiver3(hca,0,0,0,hatB0(1),hatB0(2),hatB0(3),scale,'linewidth',2)
    quiver3(hca,-scale*hatE0(1),-scale*hatE0(2),-scale*hatE0(3),hatE0(1),hatE0(2),hatE0(3),2*scale,'linewidth',2)
    
    % Label vectors
    scale = 1.7;
    text(scale*hatVExB(1),scale*hatVExB(2),scale*hatVExB(3),'v_{ExB}','fontsize',14)
    text(scale*hatVi(1),scale*hatVi(2),scale*hatVi(3),'v_{i}','fontsize',14)
    text(scale*hatVe(1),scale*hatVe(2),scale*hatVe(3),'v_{e}','fontsize',14)
    text(scale*hatB0(1),scale*hatB0(2),scale*hatB0(3),'B_{0}','fontsize',14)
    text(scale*hatE0(1),scale*hatE0(2),scale*hatE0(3),'E_{0}','fontsize',14)
    hold(hca,'off');
    hold off
end
% Set  viewing angles
for ii = 9:12; view(h(ii),hatB0); end
for ii = 5:8; view(h(ii),hatB0); end



