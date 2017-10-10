% Load electron distributions
% see mms_2015Oct16.plot_particle_distribution_1 for examples


%tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:33:50.00Z');

% Load B field
c_eval('Bvec? = dmpaB1/dmpaB1.abs;',ic);
c_eval('Bvec = Bvec?;',ic);

% ExB velocity in eV
mp = 1.6726e-27;
e = 1.6022e-19;
c_eval('vExBeV? = vExB?.abs*vExB?.abs*1e6*0.5*mp/e;',ic);
%c_eval('vExBeV? = mms.fftbandpass(vExBeV?,0,0.5);',ic);

% Ion velocity in eV
c_eval('vieV? = vi?.abs*vi?.abs*1e6*0.5*mp/e;',ic);

% Electron velocity in eV
c_eval('veeV? = ve?brst.abs*ve?brst.abs*1e6*0.5*mp/e;',ic);


% Electron moments
if 0
tmpDataObj = dataobj('data/mms1_fpi_brst_l1b_des-moms_20151001065000_v0.2.0.cdf');
ne = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_numberDensity'));
dtmoments = ne.time(2)-ne.time(1);
fsmoments = 1/dtmoments;
ne = irf_filt(ne,0,fsmoments/4,fsmoments,3);
TeX = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_TempXX'));
TeY = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_TempYY'));
TeZ = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_TempZZ'));
Te = TSeries(TeX.time,[TeX.data TeY.data TeZ.data],'to',1);
Te = irf_filt(Te,0,fsmoments/4,fsmoments,3);
end

% Electron distributions
c_eval('edist = desDist?;',ic);
c_eval('idist = disDist?;',ic);

% Define angles
dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);

% Define new PAD arrays
PSDomni = zeros(length(edist.time),length(energy));
PSDpar = PSDomni;
PSDperp = PSDomni;
PSDapar = PSDomni;
PSDpartemp = zeros(length(energy),length(phi),length(theta));
PSDperptemp = zeros(length(energy),length(phi),length(theta));
PSDapartemp = zeros(length(energy),length(phi),length(theta));

% Electron analysis
for ii = 1:length(edist.time);
    disttemp = squeeze(edist.data(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
    [~,tB] = min(abs(Bvec.time-edist.time(ii)));
    Bvecs = Bvec.data(tB,:);
    thetab = acosd(x*Bvecs(1)+y*Bvecs(2)+z*Bvecs(3));
    pospar = ones(length(phi),length(theta)); 
    pospar(find(thetab > 15)) = NaN;
    posperp = ones(length(phi),length(theta)); 
    posperp(find(thetab < 82.5)) = NaN;
    posperp(find(thetab > 97.5)) = NaN;
    posapar = ones(length(phi),length(theta)); 
    posapar(find(thetab < 165)) = NaN;    
    for kk = 1:length(energy);
        PSDpartemp(kk,:,:)  = squeeze(disttemp(kk,:,:)).*pospar;
        PSDperptemp(kk,:,:) = squeeze(disttemp(kk,:,:)).*posperp;
        PSDapartemp(kk,:,:) = squeeze(disttemp(kk,:,:)).*posapar;
    end
    PSDpar(ii,:) =  squeeze(irf.nanmean(irf.nanmean(PSDpartemp,3),2));
    PSDperp(ii,:) = squeeze(irf.nanmean(irf.nanmean(PSDperptemp,3),2));
    PSDapar(ii,:) = squeeze(irf.nanmean(irf.nanmean(PSDapartemp,3),2));
end

% Ion analysis
PSDiomni = zeros(length(idist.time),length(energy));
for ii = 1:length(idist.time);
    disttemp = squeeze(idist.data(ii,:,:,:));
    PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
end

energyspec = ones(length(edist.time),1)*energy;
efluxomni = PSDomni.*energyspec.^2;

energyspec = ones(length(idist.time),1)*energy;
ifluxomni = PSDiomni.*energyspec.^2;

specomni=struct('t',edist.time.epochUnix);
specomni.p = double(PSDomni)*1e30;
specomni.p_label={'f_e','(s^{3} km^{-6})'};
specomni.f_label={''};
specomni.f = single(energy);

specfomni=struct('t',edist.time.epochUnix);
specfomni.p = double(efluxomni);
specfomni.p_label={'e flux','(au)'};
specfomni.f_label={''};
specfomni.f = single(energy);

specfiomni=struct('t',idist.time.epochUnix);
specfiomni.p = double(ifluxomni);
specfiomni.p_label={'i flux','(au)'};
specfiomni.f_label={''};
specfiomni.f = single(energy);

specpar =struct('t',edist.time.epochUnix);
specpar.p = double(PSDpar);
specpar.p_label={'f_e par','(s^{3} m^{-6})'};
specpar.f_label={''};
specpar.f = single(energy);

specapar =struct('t',edist.time.epochUnix);
specapar.p = double(PSDapar);
specapar.p_label={'f_e perp','(s^{3} m^{-6})'};
specapar.f = single(energy);

specperp =struct('t',edist.time.epochUnix);
specperp.p = double(PSDperp);
specperp.p_label={'f_e apar','(s^{3} m^{-6})'};
specperp.f_label={''};
specperp.f = single(energy);

PSDparapar = PSDpar./PSDapar;
specparapar =struct('t',edist.time.epochUnix);
specparapar.p = double(PSDparapar);
specparapar.p_label={'f_{par}/f_{apar}'};
specparapar.f_label={''};
specparapar.f = single(energy);

PSDparperp = (PSDpar+PSDapar)./(2*PSDperp);
specparperp =struct('t',edist.time.epochUnix);
specparperp.p = double(PSDparperp);
specparperp.p_label={'(f_{par}+f_{apar})/(2 f_{perp})'};
specparperp.f_label={''};
specparperp.f = single(energy);

%%
h=irf_plot(7,'newfigure');

c_eval('dmpaB = dmpaB?.tlim(tint);',ic)
hca=irf_panel('B');
irf_plot(hca,dmpaB);
ylabel(hca,'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca=irf_panel('vi');
c_eval('irf_plot(hca,vieV?);',ic)
ylabel(hca,'v_{i} (km s^{-1})','Interpreter','tex');
irf_legend(hca,{'V_{x}','V_{y}','V_{z}'},[0.5 0.1])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca=irf_panel('n');
c_eval('irf_plot(hca,ne?brst);',ic)
hold(hca,'on');
c_eval('irf_plot(hca,ni?brst,''b'');',ic)
hold(h(3),'off');
set(hca,'yscale','log');
irf_zoom(hca,'y',[1 300]);
ylabel(hca,'n (cm^{-3})','Interpreter','tex');
irf_legend(hca,{'n_{e}','n_{i}'},[0.5 0.9]);
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca=irf_panel('idist');
irf_spectrogram(hca,specfiomni,'log','donotfitcolorbarlabel');
hold(hca,'on');
c_eval('irf_plot(hca,vieV?);',ic);
c_eval('irf_plot(hca,vExBeV?,''w'');',ic);
hold(hca,'off');
irf_legend(hca,'(d)',[0.99 0.98],'color','w','fontsize',12);
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
irf_legend(hca,{'v_i'},[0.5 0.9]);
irf_legend(hca,'v_{ExB}',[0.6 0.9],'color','w');
ylabel(hca,'E_{i} (eV)','fontsize',12,'Interpreter','tex');

hca=irf_panel('edist');
irf_spectrogram(hca,specfomni,'log','donotfitcolorbarlabel');
hold(hca,'on');
c_eval('irf_plot(hca,Te?brst)',ic)
hold(hca,'off');
irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.5 0.1])
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%caxis(hca,[-15 -5])
ylabel(hca,'E_{e} (eV)','fontsize',12,'Interpreter','tex');

hca=irf_panel('edistparperp');
irf_spectrogram(hca,specparperp,'log','donotfitcolorbarlabel');
irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,[-2 2])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,'E (eV)','fontsize',12);
colormap(hca,cn.cmap('bluered3'))

hca=irf_panel('edistparapar');
irf_spectrogram(hca,specparapar,'log','donotfitcolorbarlabel');
irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,[-2 2])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,'E (eV)','fontsize',12);
colormap(hca,cn.cmap('bluered3'))


title(h(1),'MMS 1')

%tints = irf.tint('2015-10-01T06:53:30.00Z/2015-10-01T06:54:10.00Z');
irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',tint);
set(h(1:7),'fontsize',12);

%for ii = 1: numel(h);
%    h(ii).Position(3) = h(ii).Position(3)*0.85;
%    grid(h(ii),'off')
%end

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','overviewedists1s.png');