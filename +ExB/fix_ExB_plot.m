% Load all variables
loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/';
matToLoad = '20140509T030211-2007-08-31-200V_lr0_lz5-2898280.mat';   
matToLoad = '20140509T060153-2007-08-31-200V_lr0_lz5-1448902.mat';
matToLoad = '20140517T080958-2007-08-31-200V_lr0_lz5_1300eV-3260494';
%matToLoad = '20140413T085457-2007-08-31_hotterpar2-3895000.mat';
%matToLoad = '20140409T222215_Tao2011_3051082.mat';
varToLoad = {'rSurf','zSurf','phi0','lr','lz','B0','nr','nz','rlim',...
             'zlim','meanVrz','zGrid','rGrid','Tper'};

for kk = 1:numel(varToLoad)
    %load([loadPath toLoad],varToLoad{kk});
    %if ~exist(varToLoad{kk},'var'); 
        load([loadPath,matToLoad],varToLoad{kk}); 
    %end
end

% Electric field and ExB-drift
% rz
[RS,ZS] = cn.meshgrid(rSurf,zSurf);
Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
ExBdrift = - Er/(B0*1e-9)*1e-3; % km/s
vlim = max(max(abs(ExBdrift)));

% Subtract noise level from simulated ExB-drift
% Take noise to be the mean between 10 and 20 and 80 and 90%.
limits = [0.05 0.25 0.75 0.95];
zlimits = -zlim+zlim*limits*2;
int1 = find(zSurf>(zlimits(1)) & zSurf<(zlimits(2)));
int2 = find(zSurf>(zlimits(3)) & zSurf<(zlimits(4)));
Vaz_noise = repmat(mean(meanVrz(:,[int1, int2]),2),1,nz);
no_noise = meanVrz-Vaz_noise;

% Set up figure
fig = figure(28); 
set(fig,'position',[1 1 838 955]);
nPanels = 5;
for k = 1:nPanels; h(k) = subplot(nPanels,1,k); end
isub = 1;
irf_colormap('poynting');

if 1 % Illustrate the expected ExB drift           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);
    surf(hca,rGrid,zGrid,surfa',ExBdrift')
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')    
    caxis(hca,2*vlim*[-1 1]);
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]);        
    title(hca,['ExB-drift, file: ' parseunderline(matToLoad)])
    box(hca,'on')
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);    
    surf(hca,rGrid,zGrid,surfa',meanVrz'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'km/s');    
    caxis(hca,2*vlim*[-1 1])           
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,'Average azimuthal velocity')
    box(hca,'on')    
    % Add lines where the mean is taken from within
    hold(hca,'on')
    plot(hca,[0 rlim],[1 1]*zSurf(int1(1)),'k')
    plot(hca,[0 rlim],[1 1]*zSurf(int1(end)),'k')
    plot(hca,[0 rlim],[1 1]*zSurf(int2(1)),'k')
    plot(hca,[0 rlim],[1 1]*zSurf(int2(end)),'k')
    hold(hca,'off')
end
if 1 % Plot average azimuthal velocity with subtracted noise
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);    
    surf(hca,rGrid,zGrid,surfa',Vaz_noise'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'km/s');    
    caxis(hca,2*vlim*[-1 1])           
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,'The noise')
    box(hca,'on')
end
if 1 % Plot subtracted noise
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);    
    surf(hca,rGrid,zGrid,surfa',no_noise'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'km/s');    
    caxis(hca,2*vlim*[-1 1])           
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,'Average azimuthal velocity, noise subtracted')
    box(hca,'on')
end
if 1 % Plot difference between model ExB-drift and average azimuthal velocity with subtracted noise
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);    
    surf(hca,rGrid,zGrid,surfa',ExBdrift'-no_noise'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'km/s');    
    caxis(hca,2*vlim*[-1 1])           
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,{'Difference between model ExB-drift and average azimuthal velocity with subtracted noise'})
    box(hca,'on')
end

