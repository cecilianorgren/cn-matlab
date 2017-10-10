% make a 2d plot of v_ExB_real/v_ExB_theory with rc/lr and rL/lr on the two
% axes
printModel = 1;
modelnr = 21;
ExB.model;
printModel = 0;
doPlot = 0;

nGyroOrbits = 100;

nc = 25; nr = 15;
maxGyrocenter = 7;
maxGyroradius = 4;
minGyrocenter = 0.1;
minGyroradius = 0.1;

gyrocenter = linspace(minGyrocenter,maxGyrocenter,nc);
gyroradius = linspace(minGyrocenter,maxGyroradius,nr);

v_real = zeros(nc,nr);
vaz_real = zeros(nc,nr);
v_theory = zeros(nc,nr);
tic
for kk = 1:nc;  
    %fprintf([num2str(kk) '/' num2str(nc) ' '])
    fprintf([num2str(100*kk/nc,'%.0f') ' '])
    for pp = 1:nr;
        [v_real_tmp,v_theory_tmp,vaz_real_tmp] = ExB.drift_correction(modelnr,gyrocenter(kk),gyroradius(pp),nGyroOrbits,doPlot);
        v_real(kk,pp) = v_real_tmp;
        vaz_real(kk,pp) = vaz_real_tmp;
        v_theory(kk,pp) = v_theory_tmp;
    end
end
toc
t = toc;
fprintf(['Time per bin: ' num2str(t/nc/nr) '\n'])
if 0
    %% make electron distributions and see how the covarage is
    vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
    vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)
    nParticles = 10000;
    vx0 = vtper*randn(nParticles,1)/sqrt(2);
    vy0 = vtper*randn(nParticles,1)/sqrt(2);
    vper = sqrt(vx0.^2+vy0.^2);
    me = 9.10939999999e-31;
    e = 1.6022e-19;
    fce = e*B0*1e-9/me/2/pi; % Hz, take B = 10 nT as model and get the vperp from that
    rOrbit = vper/2/pi/fce;
    x0 = -maxGyrocenter+2*maxGyrocenter*rand(nParticles,1);
    y0 = -maxGyrocenter+2*maxGyrocenter*rand(nParticles,1);
    rCenter = sqrt(x0.^2+y0.^2);
    %% binning
    edgeCenter = [gyrocenter-diff(gyrocenter(1:2)) gyrocenter(end)+diff(gyrocenter(1:2))];
    edgeOrbit = [gyroradius-diff(gyroradius(1:2)) gyroradius(end)+diff(gyroradius(1:2))];
    [nnCe,binCe] = histc(rCenter,edgeCenter);
    [nnOr,binOr] = histc(rOrbit ,edgeOrbit);
    [N edges mid loc]= histcn([rOrbit rCenter],edgeOrbit,edgeCenter);

    % bins
    % for kk=1:nc
    %     [nn,bin] = histc(rCenter,[gyrocenter-diff(gyrocenter(1:2)) gyrocenter(end)+diff(gyrocenter(1:2))]);
    %     for pp=1:nr
    %         [nnOr,binOr] = histc(rOrbit(bin==pp),[gyroradius-diff(gyroradius(1:2)) gyroradius(end)+diff(gyroradius(1:2))]);
    %         bins(kk,pp) = 
    % end
end
%%
nPlots = 3;
for k = 1:nPlots; h(k) = subplot(nPlots,1,k); end
isub = 1;
doLog = 0;   
if 0    
    hca = h(isub); isub = isub + 1;
    n_ratio = N/nParticles;
    pcolor(hca,gyrocenter,gyroradius,n_ratio);
    shading(hca,'flat');
    %irf_colormap('poynting')
    cb = colorbar('peer',hca); 
    ylabel(cb,'N_{bin}/N_{tot}')    
    xlabel(hca,'gyrocenter/l_r')
    %set(hca,'clim',max(n_ratio)*[0 1])
    ylabel(hca,'gyroradius/l_r')
    title(hca,'Particle coverage')    
    cmap = irf_colormap('poynting');
    colormap(cmap)
end
if 1    
    hca = h(isub); isub = isub + 1;
    if doLog, v_ratio = sign(vaz_real./v_theory).*abs(log10(abs(vaz_real./v_theory)));
    else v_ratio = vaz_real./v_theory; end
    pcolor(hca,gyrocenter,gyroradius,v_ratio');
    shading(hca,'flat');
    %irf_colormap('poynting')
    cb = colorbar('peer',hca); 
    if doLog, ylabel(cb,'sign(v_{real}/r_{theory})*abs(log10(v_{real}/r_{theory}))')
    else ylabel(cb,'v_{gc}/v_{ExB,perfect}','fontsize',20); end
    xlabel(hca,'r_c/l_r','fontsize',20)
    set(hca,'clim',6*[-1 1])
    ylabel(hca,'\rho_c/l_r','fontsize',20)
    title(hca,'gyrocenter velocity')    
    cmap = irf_colormap('poynting');
    colormap(cmap)
end
if 1    
    hca = h(isub); isub = isub + 1;
    if doLog, v_ratio = sign(vaz_real./v_theory).*abs(log10(abs(vaz_real./v_theory)));
    else v_ratio = vaz_real./v_theory; end
    pcolor(hca,gyrocenter,gyroradius,v_real');
    shading(hca,'flat');
    %irf_colormap('poynting')
    cb = colorbar('peer',hca); 
    ylabel(cb,'v_{gc}','fontsize',20); 
    xlabel(hca,'r_c/l_r','fontsize',20)    
    ylabel(hca,'\rho_c/l_r','fontsize',20)
    title(hca,'gyrocenter velocity')    
    cmap = irf_colormap('poynting');
    colormap(cmap)
end
if 1    
    hca = h(isub); isub = isub + 1;
    if doLog, v_ratio = sign(vaz_real./v_theory).*abs(log10(abs(vaz_real./v_theory)));
    else v_ratio = vaz_real./v_theory; end
    pcolor(hca,gyrocenter,gyroradius,v_theory');
    shading(hca,'flat');
    %irf_colormap('poynting')
    cb = colorbar('peer',hca); 
    ylabel(cb,'v_{ExB,perfect}','fontsize',20); 
    xlabel(hca,'r_c/l_r','fontsize',20)    
    ylabel(hca,'\rho_c/l_r','fontsize',20)
    title(hca,'gyrocenter velocity')    
    cmap = irf_colormap('poynting');
    colormap(cmap)
end
if 0
    hca = h(isub); isub = isub + 1;
    v_ratio = v_real./v_theory;
    pcolor(hca,gyrocenter,gyroradius,v_real');
    shading(hca,'flat');
    %irf_colormap('poynting')
    cb = colorbar('peer',hca); 
    ylabel(cb,'gyrocenter velocity:   v_{real}/r_{theory}')
    xlabel(hca,'gyrocenter/l_r')
    set(hca,'clim',2*[-1 1])
    ylabel(hca,'gyroradius/l_r')
    title('v diff')
end
    



        

