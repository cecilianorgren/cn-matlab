
iclist = 1:4;
for ic = iclist
    %%
  c_eval('Bvec? = dmpaB?/dmpaB?.abs;',ic);
  c_eval('Bvec = Bvec?;',ic);


  % Electron distributions
  c_eval('eDist = desDist?;',ic);
  c_eval('iDist = disDist?;',ic);
  c_eval('eStepTable = estepTable?;',ic);
  c_eval('iStepTable = istepTable?;',ic);
  c_eval('iEnergy! = ienergy!?.data;',ic,[0 1]);
  c_eval('eEnergy! = eenergy!?.data;',ic,[0 1]);  
  c_eval('eTheta = etheta?.data;',ic,[0 1]);
  c_eval('iTheta = itheta?.data;',ic,[0 1]);
  c_eval('ePhi = ephi?.data;',ic,[0 1]);
  c_eval('iPhi = iphi?.data;',ic,[0 1]);
  
  % Rebin data to 64 energies
  c_eval('[eDist64,ePhi64,eEnergy64] = mms.psd_rebin(eDist,ephi?,eEnergy0,eEnergy1,eStepTable);',ic)
  c_eval('[iDist64,iPhi64,iEnergy64] = mms.psd_rebin(iDist,iphi?,iEnergy0,iEnergy1,iStepTable);',ic)
  
  
  eDist.data = eDist.data*1e30; % Unit conversion
  iDist.data = iDist.data*1e30; 
  eDist64.data = eDist64.data*1e30; % Unit conversion
  iDist64.data = iDist64.data*1e30; 

  energyspec64 = ones(length(eDist64.time),1)*eEnergy64;
  ienergyspec64 = ones(length(iDist64.time),1)*iEnergy64;
  
  energyspec = ones(length(eDist.time),1)*eEnergy0;
  for ii = 1:length(eDist64.time);
      if eStepTable.data(ii),
          energyspec(ii,:) = eEnergy1;
      end
  end

  energyspeci = ones(length(iDist.time),1)*iEnergy0;
  for ii = 1:length(iDist64.time);
      if iStepTable.data(ii),
          energyspeci(ii,:) = iEnergy1;
      end
  end
  
  dangle = pi/16;
  lengthphi = 32;

  z2 = ones(lengthphi,1)*sind(eTheta);
  solida = dangle*dangle*z2;
  allsolidi = zeros(size(iDist.data));
  allsolide = zeros(size(eDist.data));
  allsolidi64 = zeros(size(iDist64.data));
  allsolide64 = zeros(size(eDist64.data));

  for ii = 1:length(iDist.time);
      for jj=1:length(iEnergy0);
          allsolidi(ii,jj,:,:) = solida;
      end
  end

  for ii = 1:length(eDist.time);
      for jj=1:length(eEnergy0);
          allsolide(ii,jj,:,:) = solida;
      end
  end

  for ii = 1:length(iDist64.time);
      for jj=1:length(iEnergy64);
          allsolidi64(ii,jj,:,:) = solida;
      end
  end

  for ii = 1:length(eDist64.time);
      for jj=1:length(eEnergy64);
          allsolide64(ii,jj,:,:) = solida;
      end
  end
  distis = iDist.data.*allsolidi;
  distes = eDist.data.*allsolide;
  distis64 = iDist64.data.*allsolidi64;
  distes64 = eDist64.data.*allsolide64;

  
  % Electron analysis - OMNI
  ePSDomni = zeros(length(eDist.time),length(eEnergy0));
  for ii = 1:length(eDist.time);
      disttemp = squeeze(distes(ii,:,:,:));
      ePSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
  end

  % Ion analysis - OMNI
  iPSDomni = zeros(length(iDist.time),length(iEnergy0));
  for ii = 1:length(iDist.time);
      disttemp = squeeze(distis(ii,:,:,:));
      iPSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
  end
  
    % Electron analysis - OMNI
  ePSDomni64 = zeros(length(eDist64.time),length(eEnergy64));
  for ii = 1:length(eDist64.time);
      disttemp = squeeze(distes64(ii,:,:,:));
      ePSDomni64(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
  end

  % Ion analysis - OMNI
  iPSDomni64 = zeros(length(iDist64.time),length(iEnergy64));
  for ii = 1:length(iDist64.time);
      disttemp = squeeze(distis64(ii,:,:,:));
      iPSDomni64(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
  end

  energyspec = ones(length(eDist.time),1)*eEnergy0;
  efluxomni = ePSDomni.*energyspec.^2;
  efluxomni = efluxomni/1e6/(5.486e-4)^2/0.53707; %convert to normal units

  energyspec = ones(length(iDist.time),1)*iEnergy0;
  ifluxomni = iPSDomni.*energyspec.^2;
  ifluxomni = ifluxomni/1e6/0.53707; %convert to normal units

  energyspec64 = ones(length(eDist64.time),1)*eEnergy64;
  efluxomni64 = ePSDomni64.*energyspec64.^2;
  efluxomni64 = efluxomni64/1e6/(5.486e-4)^2/0.53707; %convert to normal units

  energyspec64 = ones(length(iDist64.time),1)*iEnergy64;
  ifluxomni64 = iPSDomni64.*energyspec64.^2;
  ifluxomni64 = ifluxomni64/1e6/0.53707; %convert to normal units
  
  eDEFomni=struct('t',eDist.time.epochUnix);
  eDEFomni.p = double(efluxomni);
  eDEFomni.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  eDEFomni.f_label={''};
  eDEFomni.f = single(eEnergy0);
  c_eval('eDEFomni? = eDEFomni;',ic)

  iDEFomni=struct('t',iDist.time.epochUnix);
  iDEFomni.p = double(ifluxomni);
  iDEFomni.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  iDEFomni.f_label={''};
  iDEFomni.f = single(iEnergy0);
  c_eval('iDEFomni? = iDEFomni;',ic)
  
  eDEFomni64=struct('t',eDist64.time.epochUnix);
  eDEFomni64.p = double(efluxomni64);
  eDEFomni64.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  eDEFomni64.f_label={''};
  eDEFomni64.f = single(eEnergy64);
  c_eval('eDEFomni64_? = eDEFomni64;',ic)

  iDEFomni64=struct('t',iDist64.time.epochUnix);
  iDEFomni64.p = double(ifluxomni64);
  iDEFomni64.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  iDEFomni64.f_label={''};
  iDEFomni64.f = single(iEnergy64);
  c_eval('iDEFomni64_? = iDEFomni64;',ic)
  
  if 0
  
  

  energyspec = ones(length(eDist.time(1:2:end)),1)*energy;
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

  specomni=struct('t',eDist.time.epochUnix);
  specomni.p = double(PSDomni)*1e30;
  specomni.p_label={'f_e','(s^{3} km^{-6})'};
  specomni.f_label={''};
  specomni.f = single(energy);
  c_eval('ePSDomni? = specomni;',ic)

  specfomni=struct('t',eDist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxomni);
  specfomni.p_label={'log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFomni? = specfomni;',ic)

  specfomni=struct('t',eDist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxpar);
  specfomni.p_label={'par log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFpar? = specfomni;',ic)

  specfomni=struct('t',eDist.time(1:2:end).epochUnix);
  specfomni.p = double(efluxapar);
  specfomni.p_label={'apar log(dEF)','keV/(cm^2 s sr keV)'};
  specfomni.f_label={'E_e (eV)'};
  specfomni.f = single(energy);
  c_eval('eDEFapar? = specfomni;',ic)

  specfomni=struct('t',eDist.time(1:2:end).epochUnix);
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

  specpar =struct('t',eDist.time.epochUnix);
  specpar.p = double(PSDpar*1e30);
  specpar.p_label={'f_e par','(s^{3} km^{-6})'};
  specpar.f_label={'E (eV)'};
  specpar.f = single(energy);
  c_eval('ePSDpar? = specpar;',ic)

  specapar =struct('t',eDist.time.epochUnix);
  specapar.p = double(PSDapar*1e30);
  specapar.p_label={'f_e perp','(s^{3} km^{-6})'};
  specapar.f = single(energy);
  c_eval('ePSDapar? = specapar;',ic)

  specperp =struct('t',eDist.time.epochUnix);
  specperp.p = double(PSDperp*1e30);
  specperp.p_label={'f_e perp','(s^{3} km^{-6})'};
  specperp.f_label={'E (eV)'};
  specperp.f = single(energy);
  c_eval('ePSDperp? = specperp;',ic)

  PSDparapar = PSDpar./PSDapar;
  specparapar =struct('t',eDist.time.epochUnix);
  specparapar.p = double(PSDparapar);
  specparapar.p_label={'f_{par}/f_{apar}'};
  specparapar.f_label={''};
  specparapar.f = single(energy);
  c_eval('ePSDparapar? = specparapar;',ic)

  PSDparperp = (PSDpar+PSDapar)./(2*PSDperp);
  specparperp =struct('t',eDist.time.epochUnix);
  specparperp.p = double(PSDparperp);
  specparperp.p_label={'(f_{par}+f_{apar})/(2 f_{perp})'};
  specparperp.f_label={''};
  specparperp.f = single(energy);
  c_eval('ePSDparperp? = specparperp;',ic)
  end
end
