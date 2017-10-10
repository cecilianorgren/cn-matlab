% Load electron distributions
if 1,
tmpDataObj = dataobj('data/mms1_fpi_brst_l1b_des-dist_20151001065000_v0.2.0.cdf');
diste = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_brstSkyMap_dist'));
tmpDataObj = dataobj('data/mms1_fpi_brst_l1b_dis-dist_20151001065000_v0.2.0.cdf');
disti = mms.variable2ts(get_variable(tmpDataObj,'mms1_dis_brstSkyMap_dist'));
end

tint = irf.tint(diste.time.start.utc,diste.time.stop.utc);

% Load B field
tmpDataObjB = dataobj('data/mms1_dfg_srvy_ql_20151001_v1.2.1.cdf');
Bxyz = mms.variable2ts(get_variable(tmpDataObjB,'mms1_dfg_srvy_dmpa'));
Bxyz = Bxyz.tlim(tint);
Bmag = Bxyz.abs.data;
Bvec = Bxyz.data./([Bmag Bmag Bmag]);
Bvec = TSeries(Bxyz.time,Bvec,'to',1);

% Load Fast electric field data
tmpDataObj = dataobj('data/mms1_edp_fast_ql_dce2d_20151001000000_v0.2.2.cdf');
Exyz = mms.variable2ts(get_variable(tmpDataObj,'mms1_edp_dce_xyz_dsl'));

% Calculate ExB velocity
Exyz = Exyz.resample(Bxyz);
VExB = cross(Exyz,Bxyz);
VExB.data = VExB.data./[Bmag.^2 Bmag.^2 Bmag.^2];
VExBav = VExB.abs.data*1e6;
VExBav = 0.5*1.673e-27*VExBav.^2/1.6e-19;
VExBav = TSeries(VExB.time,VExBav,'to',0);
VExBav = mms.fftbandpass(VExBav,0,0.5);

% Load Ion Moments 
tmpDataObj = dataobj('data/mms1_fpi_brst_l1b_dis-moms_20151001065000_v0.2.0.cdf');
ni = mms.variable2ts(get_variable(tmpDataObj,'mms1_dis_numberDensity'));
dtmoments = ni.time(2)-ni.time(1);
fsmoments = 1/dtmoments;
ni = irf_filt(ni,0,fsmoments/4,fsmoments,3);
Vx = mms.variable2ts(get_variable(tmpDataObj,'mms1_dis_bulkX'));
Vy = mms.variable2ts(get_variable(tmpDataObj,'mms1_dis_bulkY'));
Vz = mms.variable2ts(get_variable(tmpDataObj,'mms1_dis_bulkZ'));
V = TSeries(Vx.time,[Vx.data Vy.data Vz.data],'to',1);
V = irf_filt(V,0,fsmoments/4,fsmoments,3);
Vav = V.abs.data*1000;
Vav = 0.5*1.673e-27*Vav.^2/1.6e-19;
Vav = TSeries(V.time,Vav,'to',0);




% Electron moments
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

% define angles
dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);

% Define new PAD arrays
PSDomni = zeros(length(diste.time),length(energy));
PSDpar = PSDomni;
PSDperp = PSDomni;
PSDapar = PSDomni;
PSDpartemp = zeros(length(energy),length(phi),length(theta));
PSDperptemp = zeros(length(energy),length(phi),length(theta));
PSDapartemp = zeros(length(energy),length(phi),length(theta));

% Electron analysis
for ii = 1:length(diste.time);
    disttemp = squeeze(diste.data(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
    [~,tB] = min(abs(Bvec.time-diste.time(ii)));
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
PSDiomni = zeros(length(disti.time),length(energy));
for ii = 1:length(disti.time);
    disttemp = squeeze(disti.data(ii,:,:,:));
    PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
end

energyspec = ones(length(diste.time),1)*energy;
efluxomni = PSDomni.*energyspec.^2;

energyspec = ones(length(disti.time),1)*energy;
ifluxomni = PSDiomni.*energyspec.^2;

specomni=struct('t',diste.time.epochUnix);
specomni.p = double(PSDomni)*1e30;
specomni.p_label={'f_e','(s^{3} km^{-6})'};
specomni.f_label={''};
specomni.f = single(energy);

specfomni=struct('t',diste.time.epochUnix);
specfomni.p = double(efluxomni);
specfomni.p_label={'e flux','(au)'};
specfomni.f_label={''};
specfomni.f = single(energy);

specfiomni=struct('t',disti.time.epochUnix);
specfiomni.p = double(ifluxomni);
specfiomni.p_label={'i flux','(au)'};
specfiomni.f_label={''};
specfiomni.f = single(energy);

specpar =struct('t',diste.time.epochUnix);
specpar.p = double(PSDpar);
specpar.p_label={'f_e par','(s^{3} m^{-6})'};
specpar.f_label={''};
specpar.f = single(energy);

specapar =struct('t',diste.time.epochUnix);
specapar.p = double(PSDapar);
specapar.p_label={'f_e perp','(s^{3} m^{-6})'};
specapar.f = single(energy);

specperp =struct('t',diste.time.epochUnix);
specperp.p = double(PSDperp);
specperp.p_label={'f_e apar','(s^{3} m^{-6})'};
specperp.f_label={''};
specperp.f = single(energy);

PSDparapar = PSDpar./PSDapar;
specparapar =struct('t',diste.time.epochUnix);
specparapar.p = double(PSDparapar);
specparapar.p_label={'f_{par}/f_{apar}'};
specparapar.f_label={''};
specparapar.f = single(energy);

PSDparperp = (PSDpar+PSDapar)./(2*PSDperp);
specparperp =struct('t',diste.time.epochUnix);
specparperp.p = double(PSDparperp);
specparperp.p_label={'(f_{par}+f_{apar})/(2 f_{perp})'};
specparperp.f_label={''};
specparperp.f = single(energy);


h=irf_plot(7,'newfigure');

Bxyz = Bxyz.tlim(tint);
h(1)=irf_panel('Bxyz');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Vxyz');
irf_plot(h(2),V);
ylabel(h(2),'V_{DSL} (km s^{-1})','Interpreter','tex');
irf_legend(h(2),{'V_{x}','V_{y}','V_{z}'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('ne');
irf_plot(h(3),ne);
hold(h(3),'on');
irf_plot(h(3),ni,'b');
hold(h(3),'off');
set(h(3),'yscale','log');
irf_zoom(h(3),'y',[1 300]);
ylabel(h(3),'n (cm^{-3})','Interpreter','tex');
irf_legend(h(3),{'n_{e}','n_{i}'},[0.5 0.9]);
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('idist');
irf_spectrogram(h(4),specfiomni,'log','donotfitcolorbarlabel');
hold(h(4),'on');
irf_plot(h(4),Vav);
irf_plot(h(4),VExBav,'w');
hold(h(4),'off');
irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12);
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
irf_legend(h(4),{'V'},[0.5 0.9]);
irf_legend(h(4),'V_{ExB}',[0.6 0.9],'color','w');
ylabel(h(4),'E_{i} (eV)','fontsize',12,'Interpreter','tex');

h(5)=irf_panel('edist');
irf_spectrogram(h(5),specfomni,'log','donotfitcolorbarlabel');
hold(h(5),'on');
irf_plot(h(5),Te)
hold(h(5),'off');
irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(h(5),{'T_{x}','T_{y}','T_{z}'},[0.5 0.1])
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
%caxis(h(5),[-15 -5])
ylabel(h(5),'E_{e} (eV)','fontsize',12,'Interpreter','tex');

h(6)=irf_panel('edistparapar');
irf_spectrogram(h(6),specparapar,'log','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
set(h(6),'yscale','log');
caxis(h(6),[-2 2])
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),'E (eV)','fontsize',12);

h(7)=irf_panel('edistparperp');
irf_spectrogram(h(7),specparperp,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)
set(h(7),'yscale','log');
caxis(h(7),[-2 2])
set(h(7),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(7),'E (eV)','fontsize',12);

title(h(1),'MMS 1')

%tints = irf.tint('2015-10-01T06:53:30.00Z/2015-10-01T06:54:10.00Z');
irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',tint);
set(h(1:7),'fontsize',12);


%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','overviewedists1s.png');