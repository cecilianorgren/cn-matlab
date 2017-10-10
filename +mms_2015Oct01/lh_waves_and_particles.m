% first run 
% mms.lh_loaddata 
% mms.lh_get_vph
tintEpoch = toepoch([2015 08 15 12 59 50; 2015 08 15 12 59 51])';
tint = EpochTT(irf_time(tintEpoch,'epoch>utc'));
%%

disp(['B0 = [' num2str(z(1,:)) '] * ' num2str(B0) ' nT'])
disp(['vph = [' num2str(direction) '] * ' num2str(velocity) ' km/s'])

v_ExB = cross(E,dcB.resample(E));
v_ExB = mean(v_ExB.tlim(tint).data);

c_eval('v_i = mean(vi?.tlim(tint).data);',sc)
v_i_hat = double(irf_norm(v_i));
c_eval('v_e = mean(ve?.tlim(tint).data);',sc)
v_e_hat = double(irf_norm(v_e));

%% Plot a skymap with vph and B0 inserted
% Set up coordinates for skymap plot
c_eval('desSkyMap = desSkyMap?;',sc);
c_eval('disSkyMap = disSkyMap?;',sc);
r = 1;
phi_edges = linspace(0,2*pi,size(desSkyMap.data,3)+1);
theta_edges = linspace(0,pi,size(desSkyMap.data,4)+1);
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = r*sin(THETA).*cos(PHI);
Y = r*sin(THETA).*sin(PHI);
Z = r*cos(THETA);


energyLevel = 10;
[tId,~] = desSkyMap.time.tlim(tint); 
C = squeeze(nansum(desSkyMap.data(tId,energyLevel,:,:),1))';
hs = surf(X,Y,Z,C);
axis square
axis equal
hc = colorbar;
hca=gca;
hca.XLabel.String = 'X';
hca.YLabel.String = 'Y';
hca.ZLabel.String = 'Z';
hca.Title.String = tint.utc;

hold on;
scale = 1.5;
v_ExB_hat = double(irf_norm(v_ExB));
quiver3(0,0,0,x(i_dir,1),x(i_dir,2),x(i_dir,3),scale,'linewidth',2)
quiver3(0,0,0,y(i_dir,1),y(i_dir,2),y(i_dir,3),scale,'linewidth',2)
quiver3(0,0,0,z(i_dir,1),z(i_dir,2),z(i_dir,3),scale,'linewidth',2)
quiver3(0,0,0,v_ExB_hat(1),v_ExB_hat(2),v_ExB_hat(3),scale,'linewidth',2)
quiver3(0,0,0,v_i_hat(1),v_i_hat(2),v_i_hat(3),scale,'linewidth',2)

scale = 1.7;
text(scale*x(i_dir,1),scale*x(i_dir,2),scale*x(i_dir,3),'v_{ph}','fontsize',14)
text(scale*y(i_dir,1),scale*y(i_dir,2),scale*y(i_dir,3),'n','fontsize',14)
text(scale*z(i_dir,1),scale*z(i_dir,2),scale*z(i_dir,3),'B','fontsize',14)
text(scale*v_ExB_hat(1),scale*v_ExB_hat(2),scale*v_ExB_hat(3),'v_{ExB}','fontsize',14)

%% Plot a skymap for one time and several energy levels with vph and B0 inserted
% Set up coordinates for skymap plot
c_eval('desSkyMap = desSkyMap?;',sc);
c_eval('disSkyMap = disSkyMap?;',sc);
skyMap = disSkyMap;
r = 1;
phi_edges = linspace(0,2*pi,size(skyMap.data,3)+1);
theta_edges = linspace(0,pi,size(skyMap.data,4)+1);
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = -r*sin(THETA).*cos(PHI);
Y = -r*sin(THETA).*sin(PHI);
Z = -r*cos(THETA);

% Plot a few energy levels on a skymap
nrows = 3;
ncols = 4;
nEnergyLevels = size(skyMap.data,2);
maxEnergyLevel = 32;
levelStep = ceil((maxEnergyLevel)/nrows/ncols);
nPlots = floor(maxEnergyLevel/levelStep);

[tId,~] = skyMap.time.tlim(tint); 
energyLevels = 1:levelStep:maxEnergyLevel; 
%energyLevels = [15:1:26];
nPlots = numel(energyLevels);
for k = 1:nPlots
    energyLevel = energyLevels(k);%1+(k-1)*levelStep;
    hca = subplot(nrows,ncols,k);
    h(k) = hca;
    C = squeeze(nanmean(skyMap.data(tId,energyLevel,:,:),1))';
    hs = surf(hca,X,Y,Z,C);
    hc = colorbar('peer',hca);
    axis square
    axis equal

    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(energyLevel)]};
    hca.Title.String = titleString;
    %view([-1 1 1])    
    shading flat;
    hold(hca,'on');
    hold on
    scale = 1.5;
    v_ExB_hat = double(irf_norm(v_ExB));
    quiver3(0,0,0,x(i_dir,1),x(i_dir,2),x(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,y(i_dir,1),y(i_dir,2),y(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,z(i_dir,1),z(i_dir,2),z(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,v_ExB_hat(1),v_ExB_hat(2),v_ExB_hat(3),scale,'linewidth',2)
    quiver3(0,0,0,v_i_hat(1),v_i_hat(2),v_i_hat(3),scale,'linewidth',2)
    quiver3(0,0,0,v_e_hat(1),v_e_hat(2),v_e_hat(3),scale,'linewidth',2)
    
    scale = 1.7;
    text(scale*x(i_dir,1),scale*x(i_dir,2),scale*x(i_dir,3),'v_{ph}','fontsize',14)
    text(scale*y(i_dir,1),scale*y(i_dir,2),scale*y(i_dir,3),'n','fontsize',14)
    text(scale*z(i_dir,1),scale*z(i_dir,2),scale*z(i_dir,3),'B','fontsize',14)
    text(scale*v_ExB_hat(1),scale*v_ExB_hat(2),scale*v_ExB_hat(3),'v_{ExB}','fontsize',14)
    text(scale*v_i_hat(1),scale*v_i_hat(2),scale*v_i_hat(3),'v_{i}','fontsize',14)
    text(scale*v_i_hat(1),scale*v_e_hat(2),scale*v_e_hat(3),'v_{e}','fontsize',14)
    hold(hca,'off');
    hold off
end
for ii = 1:nPlots; view(h(ii),x(i_dir,:)+y(i_dir,:)+z(i_dir,:)); end

%% Plot a skymap for a few times with two energy levels with vph and B0 inserted
% 3 rows, electric field/potential in one
% skymap for one energy in 2nd row
% skymap for another energy in 3d row

% Set up coordinates for skymap plot
c_eval('desSkyMap = desSkyMap?;',sc);
c_eval('disSkyMap = disSkyMap?;',sc);
skyMap = disSkyMap;
r = 1;
phi_edges = linspace(0,2*pi,size(skyMap.data,3)+1);
theta_edges = linspace(0,pi,size(skyMap.data,4)+1);
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = -r*sin(THETA).*cos(PHI);
Y = -r*sin(THETA).*sin(PHI);
Z = r*cos(THETA);

% Plot a few energy levels on a skymap
nrows = 3;
ncols = 4;
nPlots = ncols;

% Choose energy levels
energyLEvel1 = 10;
energyLEvel2 = 15;

% Define time intervals


[tIdWaves,~] = skyMap.time.tlim(tint); % time for when the potential is taken

hca = subplot(nrows,ncols,1:ncols);
irf_plot(hca,{phi_B,phi_E(:,[1 1+i_v])},'comp')
irf_legend(hca,{'\phi_B','\phi_E'},[0.98 0.95])

for k = 1:nPlots    
    % Energy level 1
    hca = subplot(nrows,ncols,ncols+k); % 
    h(k) = hca;
    C = squeeze(nanmean(skyMap.data(tId,energyLevel,:,:),1))';
    hs = surf(hca,X,Y,Z,C);
    hc = colorbar('peer',hca);
    axis square
    axis equal

    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(energyLevel)]};
    hca.Title.String = titleString;
    %view([-1 1 1])    
    shading flat;
    hold(hca,'on');
    hold on
    scale = 1.5;
    v_ExB_hat = double(irf_norm(v_ExB));
    quiver3(0,0,0,x(i_dir,1),x(i_dir,2),x(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,y(i_dir,1),y(i_dir,2),y(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,z(i_dir,1),z(i_dir,2),z(i_dir,3),scale,'linewidth',2)
    quiver3(0,0,0,v_ExB_hat(1),v_ExB_hat(2),v_ExB_hat(3),scale,'linewidth',2)
    quiver3(0,0,0,v_i_hat(1),v_i_hat(2),v_i_hat(3),scale,'linewidth',2)
    quiver3(0,0,0,v_e_hat(1),v_e_hat(2),v_e_hat(3),scale,'linewidth',2)

    scale = 1.7;
    text(scale*x(i_dir,1),scale*x(i_dir,2),scale*x(i_dir,3),'v_{ph}','fontsize',14)
    text(scale*y(i_dir,1),scale*y(i_dir,2),scale*y(i_dir,3),'n','fontsize',14)
    text(scale*z(i_dir,1),scale*z(i_dir,2),scale*z(i_dir,3),'B','fontsize',14)
    text(scale*v_ExB_hat(1),scale*v_ExB_hat(2),scale*v_ExB_hat(3),'v_{ExB}','fontsize',14)
    text(scale*v_i_hat(1),scale*v_i_hat(2),scale*v_i_hat(3),'v_{i}','fontsize',14)
    text(scale*v_i_hat(1),scale*v_e_hat(2),scale*v_e_hat(3),'v_{e}','fontsize',14)
    hold(hca,'off');
    hold off
end
for ii = 1:nPlots; view(h(ii),x(i_dir,:)+y(i_dir,:)+z(i_dir,:)); end