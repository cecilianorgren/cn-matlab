%%
  iclist = 1:4;
  for ic = iclist
    %%
  c_eval('Bvec? = dmpaB?/dmpaB?.abs;',ic);
  c_eval('Bvec = Bvec?;',ic);


  % Electron distributions
  c_eval('edist = desDist?;',ic);
  c_eval('idist = disDist?;',ic);

  if 0 % Needed or not?
    edist.data = edist.data*1e30; % Unit conversion
    idist.data = idist.data*1e30;
  end

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
  %%
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
  %%
  % Ion analysis
  PSDiomni = zeros(length(idist.time),length(energy));
  for ii = 1:length(idist.time);
      disttemp = squeeze(idist.data(ii,:,:,:));
      PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3));
  end

  energyspec = ones(length(edist.time(1:2:end)),1)*energy;
  efluxomni = PSDomni(1:2:end,:).*energyspec.^2;
  efluxomni = efluxomni/1e6/(5.486e-4)^2/0.53707; %convert to normal units
  efluxomni = efluxomni*1e30; % this seems appropriate /Cecilia

  efluxpar = PSDpar(1:2:end,:).*energyspec.^2;
  efluxpar = efluxpar/1e6/(5.486e-4)^2/0.53707; %convert to normal units
  efluxpar = efluxpar*1e30; % this seems appropriate /Cecilia

  efluxapar = PSDapar(1:2:end,:).*energyspec.^2;
  efluxapar = efluxapar/1e6/(5.486e-4)^2/0.53707; %convert to normal units
  efluxapar = efluxapar*1e30; % this seems appropriate /Cecilia

  efluxperp = PSDperp(1:2:end,:).*energyspec.^2;
  efluxperp = efluxperp/1e6/(5.486e-4)^2/0.53707; %convert to normal units
  efluxperp = efluxperp*1e30; % this seems appropriate /Cecilia


  energyspec = ones(length(idist.time(1:2:end)),1)*energy;
  ifluxomni = PSDiomni(1:2:end,:).*energyspec.^2;
  ifluxomni = ifluxomni/1e6/0.53707; %convert to normal units
  ifluxomni = ifluxomni*1e30; % this seems appropriate /Cecilia

  specomni=struct('t',edist.time.epochUnix);
  specomni.p = double(PSDomni)*1e30;
  specomni.p_label={'f_e','(s^{3} km^{-6})'};
  specomni.f_label={''};
  specomni.f = single(energy);
  c_eval('ePSDomni? = specomni;',ic)

  specfomni=struct('t',edist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxomni);
  specfomni.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFomni? = specfomni;',ic)

  specfomni=struct('t',edist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxpar);
  specfomni.p_label={'par log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFpar? = specfomni;',ic)

  specfomni=struct('t',edist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxapar);
  specfomni.p_label={'apar log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFapar? = specfomni;',ic)

  specfomni=struct('t',edist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxperp);
  specfomni.p_label={'perp log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFperp? = specfomni;',ic)

  specfiomni=struct('t',idist.time(1:2:end).epochUnix);
  specfiomni.p = double(ifluxomni);
  specfiomni.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  specfiomni.f_label={'E_i (eV)'};
  specfiomni.f = single(energy);
  c_eval('iDEFomni? = specfiomni;',ic)

  specpar =struct('t',edist.time.epochUnix);
  specpar.p = double(PSDpar*1e30);
  specpar.p_label={'f_e par','(s^{3} km^{-6})'};
  specpar.f_label={'E (eV)'};
  specpar.f = single(energy);
  c_eval('ePSDpar? = specpar;',ic)

  specapar =struct('t',edist.time.epochUnix);
  specapar.p = double(PSDapar*1e30);
  specapar.p_label={'f_e perp','(s^{3} km^{-6})'};
  specapar.f = single(energy);
  c_eval('ePSDapar? = specapar;',ic)

  specperp =struct('t',edist.time.epochUnix);
  specperp.p = double(PSDperp*1e30);
  specperp.p_label={'f_e perp','(s^{3} km^{-6})'};
  specperp.f_label={'E (eV)'};
  specperp.f = single(energy);
  c_eval('ePSDperp? = specperp;',ic)

  PSDparapar = PSDpar./PSDapar;
  specparapar =struct('t',edist.time.epochUnix);
  specparapar.p = double(PSDparapar);
  specparapar.p_label={'f_{par}/f_{apar}'};
  specparapar.f_label={''};
  specparapar.f = single(energy);
  c_eval('ePSDparapar? = specparapar;',ic)

  PSDparperp = (PSDpar+PSDapar)./(2*PSDperp);
  specparperp =struct('t',edist.time.epochUnix);
  specparperp.p = double(PSDparperp);
  specparperp.p_label={'(f_{par}+f_{apar})/(2 f_{perp})'};
  specparperp.f_label={''};
  specparperp.f = single(energy);
  c_eval('ePSDparperp? = specparperp;',ic)

  end
