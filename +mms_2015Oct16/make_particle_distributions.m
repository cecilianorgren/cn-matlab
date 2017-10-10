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