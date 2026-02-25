ic=1;
if numel(ic)==4 & all(ic==1:4), do4sc = 1; else do4sc = 0; end
units = irf_units;

if do4sc 
  % Calculate curlometer current
  c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
  [Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
  gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
  gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
  gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
  
  % Calculate gradients
  gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km';
  gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';
  gseGradTe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';
  gseGradTi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';
  gseGradNe = c_4_grad('gseR?','ne?','grad');
  gseGradNi = c_4_grad('gseR?','ni?','grad');
  
  % Calculate curvature of magnetic field and electron velocity (inertial term)
  c_eval('B? = gseB?.resample(gseB1); B? = [B?.time.epochUnix double(B?.data)];',1:4)
  c_eval('R? = gseR?.resample(gseB1); R? = [R?.time.epochUnix R?.data];',1:4)
  [curvB,BB]=c_4_grad('B?','R?','curvature');
  curvB = irf.ts_vec_xyz(gseB1.time,curvB(:,2:4));
  BB = irf.ts_scalar(gseB1.time,BB(:,1:3));
  
  c_eval('Ve? = gseVe?.resample(gseVe1); Ve? = [Ve?.time.epochUnix double(Ve?.data)];',1:4)
  c_eval('R? = gseR?.resample(gseVe1); R? = [R?.time.epochUnix R?.data];',1:4)
  [gseCurvVe,VVE]=c_4_grad('Ve?','R?','curvature');
  gseCurvVe = irf.ts_vec_xyz(gseVe1.time,gseCurvVe(:,2:4));
  VVE = irf.ts_scalar(gseVe1.time,VVE(:,1:3));
  
  % Inertial term, not complete
  divVe = c_4_grad('gseR?','gseVe?','div'); divVe.name = 'div Ve';
end

% Calculate currents from moments and ExB
c_eval('gseJExB? = -units.e*ne?.resample(gseVExB?)*gseVExB?*1e3*1e6*1e9; gseJExB?.units = ''nA/m^2''; gseJExB?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

% Calculate convective electric fields
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

% ExB velocity
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3;',ic) % km/s

% JxB
c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?));',ic)

% Non-ideal electric field, E+VexB
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

% Decompose into parallel and perpendicular components
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

if do4sc
  [gseGradPepar,gseGradPeperp] = irf_dec_parperp(gseAvB,gseGradPe); gseGradPepar.name = 'div Pe par'; gseGradPeperp.name = 'div Pe par';
end

% Field aligned electric and magnetic fields
c_eval('facB? = irf_convert_fac(gseB?,gseB?,[1 0 0]);',ic)
c_eval('facB?scm = irf_convert_fac(gseB?scm,gseB?,[1 0 0]);',ic)
c_eval('facE? = irf_convert_fac(gseE?,gseB?,[1 0 0]);',ic)
%% Calculate some physical parameters with irf_plasma_calc
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)
 
% Speeds
c_eval('matB? = gseB?.abs.data,ne?.resample(gseB?.time).data;',ic)
c_eval('matParTe? = facTe?.xx.resample(gseB?.time).data;',ic)
c_eval('matParTi? = facTi?.xx.resample(gseB?.time).data;',ic)
c_eval('matPerTe? = (facTe?.yy.resample(gseB?.time).data + facTe?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matPerTi? = (facTi?.yy.resample(gseB?.time).data + facTi?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matTe? = facTe?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matTi? = facTi?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matNe? = ne?.resample(gseB?.time).data;',ic)

c_eval('vte?perp = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Vte''); vte?perp = irf.ts_scalar(gseB?.time,vte?perp)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte?par = irf_plasma_calc(matB?,matNe?,0,matParTe?,matTi?,''Vte''); vte?par = irf.ts_scalar(gseB?.time,vte?par)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Vte''); vte? = irf.ts_scalar(gseB?.time,vte?)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vtp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matPerTi?,''Vtp''); vtp? = irf.ts_scalar(gseB?.time,vtp?)*1e-3; vtp?.units = ''km/s'';',ic)
c_eval('vA? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Va''); vA? = irf.ts_scalar(gseB?.time,vA?)*1e-3; vA?.units = ''km/s'';',ic)

% Frequencies
c_eval('flh? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Flh''); flh? = irf.ts_scalar(gseB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Fce''); fce? = irf.ts_scalar(gseB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fcp''); fcp? = irf.ts_scalar(gseB?.time,fcp?);',ic)
c_eval('fpe? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpe''); fpe? = irf.ts_scalar(gseB?.time,fpe?);',ic)
c_eval('fpp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpp''); fpp? = irf.ts_scalar(gseB?.time,fpp?);',ic)

% Length scales
c_eval('Lp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Li''); Lp? = irf.ts_scalar(gseB?.time,Lp?)*1e-3; Lp?.units = ''km''; Lp?.name=''p inertial length'';',ic)
c_eval('Le? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Le''); Le? = irf.ts_scalar(gseB?.time,Le?)*1e-3; Le?.units = ''km''; Le?.name=''e inertial length'';',ic)
c_eval('Ld? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Ld''); Ld? = irf.ts_scalar(gseB?.time,Ld?)*1e-3; Ld?.units = ''km''; Ld?.name=''Debye length'';',ic)
c_eval('re? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Roe''); re? = irf.ts_scalar(gseB?.time,re?)*1e-3; re?.units = ''km''; re?.name=''e gyroradius'';',ic)
c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gseB?.time,rp?)*1e-3; rp?.units = ''km''; re?.name=''p gyroradius'';',ic)

c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)



% Pitchangle distribution
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)



%% Wave analysis
if 0
  %%
c_eval('wavE? = irf_wavelet(gseE?.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
%%
tintPol = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:52.00Z');
c_eval('polarization? = irf_ebsp(gseE?.tlim(tintPol),gseB?scm.tlim(tintPol),gseB?.tlim(tintPol),gseB?.tlim(tintPol),gseR?brsttime.tlim(tintPol),[10 2200],''polarization'',''fac'');',1)
end




%% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

% Compute electron Mach number
c_eval('Me?perp = gseVe?perp.abs/vte?perp;',ic);
c_eval('Me?par = gseVe?par.abs/vte?par;',ic);

% Compute current density and J.E
c_eval('EdotJ? = gseE?.resample(gseJ?).dot(gseJ?)/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?par = gseE?par.resample(gseJ?par)*gseJ?par/1000; EdotJ?par.units = ''nW/m^3''; EdotJ?par.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?perp = gseE?perp.resample(gseJ?perp).dot(gseJ?perp)/1000; EdotJ?perp.units = ''nW/m^3''; EdotJ?perp.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('RedotJ? = gseEVexB?.resample(gseJ?).dot(gseJ?)/1000; RedotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

% Calculate epsilon and delta parameters
c_eval('oce? = fce?*2*pi;',ic)
c_eval('EdotVe? = gseE?.resample(gseVe?).dot(gseVe?);',ic);
c_eval('epsilone? = abs(6*pi*EdotVe?/(oce?.resample(gseVe?)*(facTe?.trace)));',ic);


c_eval('deltae? = gseVexB?/(gseVe?perp.abs*gseB?.resample(gseVe?).abs*1e-9);',ic);
%c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

gseR0 = (gseR1.resample(gseR1.time)+gseR2.resample(gseR1.time)+gseR3.resample(gseR1.time)+gseR4.resample(gseR1.time))/4;
c_eval('gseRR? = gseR?-gseR0; gseRR? = gseRR?.resample(irf_time(''2015-11-12T07:19:21.000Z'',''utc>epochTT'')).data;',ic)


% Plasma beta and magnetic pressure
c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)


%% Average properties
if do4sc
  gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 
  avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
  gseAvE = (gseE1+gseE2.resample(gseE1.time)+gseE3.resample(gseE1.time)+gseE4.resample(gseE1.time))/4; 
  gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; 
  gseAvVeperp = (gseVe1perp+gseVe2perp.resample(gseVe1perp.time)+gseVe3perp.resample(gseVe1perp.time)+gseVe4perp.resample(gseVe1perp.time))/4; 
  gseAvVi = (gseVi1+gseVi2.resample(gseVi1.time)+gseVi3.resample(gseVi1.time)+gseVi4.resample(gseVi1.time))/4; 
  gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; 
  gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.units = gseJ1.units;
  gseAvVExB = (gseVExB1+gseVExB2.resample(gseVExB1.time)+gseVExB3.resample(gseVExB1.time)+gseVExB4.resample(gseVExB1.time))/4; 
  gseAvVexB = (gseVexB1+gseVexB2.resample(gseVexB1.time)+gseVexB3.resample(gseVexB1.time)+gseVexB4.resample(gseVexB1.time))/4; 
  gseAvVixB = (gseVixB1+gseVixB2.resample(gseVixB1.time)+gseVixB3.resample(gseVixB1.time)+gseVixB4.resample(gseVixB1.time))/4; 
  gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4;
end


% Ohm's law terms
%gseAvEVexB = (gseEVexB1 + gseEVexB2.resample(gseEVexB1) + gseEVexB3.resample(gseEVexB1) + gseEVexB4.resample(gseEVexB1))/4;
if do4sc
  gseOhmGradPe = gseGradPe/avNe.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseOhmGradPe.units = 'mV/m';
  gseOhmVexB = gseAvVexB; gseOhmVexB.units = 'mV/m';
  gseOhmVixB = gseAvVixB; gseOhmVixB.units = 'mV/m';
  gseOhmJxB_a = gseAvJ.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_a.units = 'mV/m';
  gseOhmJxB_b = gseJcurl.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_b.units = 'mV/m';
  gseOhmJxB_c = (gseJxB1/ne1+gseJxB2.resample(gseJxB1.time)/ne2+gseJxB3.resample(gseJxB1.time)/ne3+gseJxB4.resample(gseJxB1.time)/ne4)/4/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_c.units = 'mV/m'; 
  gseOhmJxB = gseOhmJxB_c;
end




if 0
%% MVA: Rotate into lmn coordinates
ic=1:4;

avPe = (gsePe1.trace+gsePe2.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1))/3/4;
LgradP = avPe/gseGradPe.abs; LgradP.name = 'L_{P}';

% Proxy for magnetic field slippage (Scudder2105)
avFce = (fce1+fce2.resample(fce1)+fce3.resample(fce1)+fce4.resample(fce1))/4;
avRhoe = (re1+re2.resample(re1)+re3.resample(re1)+re4.resample(re1))/4;
Y = avFce*(avRhoe/LgradP).^2; Y.units = 's^-1';
Ynodim = (avRhoe/LgradP); Y.units = 's^-1';

facAvTe = (facTe1+facTe2.resample(facTe1)+facTe3.resample(facTe1)+facTe4.resample(facTe1))/4;

% Curl of E+VexB (Scudder2105)
gseRotRe = mms_2015Oct16.rotRe(gseR1,gseR2,gseR3,gseR4,gseEVexB1,gseEVexB2,gseEVexB3,gseEVexB4);

% Set up coordinate system
% Default
%[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
%L = v(1,:); M = v(2,:); N = v(3,:);
%coordLabels = {'L','M','N'};
%lmn = [N;-M;L];

coordSystem = 6;
switch coordSystem % Choose another
  case 1 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
    L = v(1,:); M = -v(2,:); N = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 2 % minimum variance of gradPe
    [out,l,v] = irf_minvar(gseGradPe.tlim(irf.tint('2015-11-12T07:19:20.484Z/2015-11-12T07:19:22.034Z')));
    N = -v(1,:); L = -v(2,:); M = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 3 % B in magnetosphere
    L = gseB1.resample(irf_time('2015-11-12T07:19:05.000Z','utc>epochtt')); L = L.data/norm(L.data);
    M = cross(L,cross([0 1 0],L)); M = M/norm(M);
    N = cross(L,M);
    lmn = [L;M;N];
  case 4 % N: minimum variance of J
    [out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-11-12T07:19:20.516Z/2015-11-12T07:19:21.826Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'N','-M','L'};
    lmn = [L;M;N];        
  case 5 % N: minimum variance of B,'<Bn>=0' 
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')),'<Bn>=0');
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 6 % N: minimum variance of B,'<Bn>=0' 
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')),'<Bn>=0');
    L = v(1,:); M = -v(2,:); N = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 22 % N: minimum variance of J
    [out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-10-16T10:33:23.000Z/2015-10-16T10:33:32.000Z')));
    L = -v(2,:); M = -v(1,:); N = v(3,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
    [out,l,v] = irf_minvar(gseJcurl.tlim(irf.tint('2015-10-16T10:33:22.595Z/2015-10-16T10:33:31.284Z')));
    L = -v(1,:); M = v(2,:); N = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];        
  case 33 % N: maximum variance of E
    tint = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.400Z');
    [out,l,v] = irf_minvar(gseE3.tlim(tint));
    L = v(2,:); M = -v(3,:); N = -v(1,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
  case 44 % N: magnetosheath side normal derived from mms1 and mms4
    gseVec14 = gseR4-gseR1; gseVec14 = gseVec14.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec14.time,M);
    gseNorm14 = gseVec14.cross(gseM);
    gseNormalMSH = gseNorm14/gseNorm14.abs;

    N = -gseNormalMSH.data;
    M = M;
    L = cross(M,N);
  case 55 % N: magnetosphere side normal derived from mms3 and mms4
    gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec34.time,M);
    gseNorm34 = gseVec34.cross(gseM);
    gseNormalMSP = gseNorm34/gseNorm34.abs;

    N = -gseNormalMSP.data;
    M = M;
    L = cross(M,N);
end

% Rotate data
c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';')
c_eval('mvaB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(L).data gseB?.dot(M).data gseB?.dot(N).data]); mvaB?.name = ''B LMN'';')
c_eval('mvaB?_ = gseB?*lmn''; mvaB?.name = ''B LMN'';')
c_eval('mvaE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(L).data gseE?.dot(M).data gseE?.dot(N).data]); mvaE?.name = ''E LMN'';')
c_eval('mvaVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(L).data gseVe?.dot(M).data gseVe?.dot(N).data]); mvaVe?.name = ''Ve LMN'';')
c_eval('mvaVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(L).data gseVi?.dot(M).data gseVi?.dot(N).data]); mvaVi?.name = ''Vi LMN'';')
c_eval('mvaJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(L).data gseJ?.dot(M).data gseJ?.dot(N).data]); mvaJ?.units = gseJ?.units; mvaJ?.name = ''J LMN'';')
c_eval('mvaJe? = irf.ts_vec_xyz(gseJe?.time,[gseJe?.dot(L).data gseJe?.dot(M).data gseJe?.dot(N).data]); mvaJe?.units = gseJe?.units; mvaJe?.name = ''Je LMN'';')
c_eval('mvaJi? = irf.ts_vec_xyz(gseJi?.time,[gseJi?.dot(L).data gseJi?.dot(M).data gseJi?.dot(N).data]); mvaJi?.units = gseJi?.units; mvaJi?.name = ''Ji LMN'';')
mvaJcurl = irf.ts_vec_xyz(gseJcurl.time,[gseJcurl.dot(L).data gseJcurl.dot(M).data gseJcurl.dot(N).data]); mvaJcurl.units = gseJcurl.units;
c_eval('mvaPe? = mms.rotate_tensor(gsePe?,''rot'',L,M,N); mvaPe? = irf.ts_tensor_xyz(mvaPe?.time,mvaPe?.data); mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaPi? = mms.rotate_tensor(gsePi?,''rot'',L,M,N);',ic)
c_eval('mvaTi? = mms.rotate_tensor(gseTi?,''rot'',L,M,N);',ic)
c_eval('mvaTe? = mms.rotate_tensor(gseTe?,''rot'',L,M,N); mvaTe? = irf.ts_tensor_xyz(mvaTe?.time,mvaTe?.data); mvaTe?.units = gseTe?.units;',ic)
c_eval('mvaVexB? = irf.ts_vec_xyz(gseVexB?.time,[gseVexB?.dot(L).data gseVexB?.dot(M).data gseVexB?.dot(N).data]); mvaVexB?.units = ''mV/m'';')
c_eval('mvaVixB? =  irf.ts_vec_xyz(gseVixB?.time,[gseVixB?.dot(L).data gseVixB?.dot(M).data gseVixB?.dot(N).data]); mvaVixB?.units = ''mV/m'';')
c_eval('mvaEVexB? =  irf.ts_vec_xyz(gseEVexB?.time,[gseEVexB?.dot(L).data gseEVexB?.dot(M).data gseEVexB?.dot(N).data]); mvaEVexB?.units = ''mV/m'';')
%mvaVDe =  irf.ts_vec_xyz(vDe.time,[vDe.dot(L).data vDe.dot(M).data vDe.dot(N).data]); mvaVDe.units = '';
mvaAvJ =  irf.ts_vec_xyz(gseAvJ.time,[gseAvJ.dot(L).data gseAvJ.dot(M).data gseAvJ.dot(N).data]); mvaAvJ.units = 'nA/m^2';
c_eval('mvaJxB? = mvaJ?.cross(mvaB?.resample(mvaJ?.time));')
c_eval('mvaVExB? =  gseVExB?*lmn'';')
c_eval('mvaVe?par = gseVe?par;')
c_eval('mvaVe?perp = gseVe?perp*lmn''; mvaVe?perp.name = ''Ve perp lmn'';')
c_eval('mvaE?par = gseE?par;')
c_eval('mvaE?perp = gseE?perp*lmn''; mvaE?perp.name = ''E perp lmn'';')
%c_eval('mvaE?fastpar = gseE?fastpar;')
%c_eval('mvaE?fastperp = irf.ts_vec_xyz(gseE?fastperp.time,[gseE?fastperp.dot(L).data gseE?fastperp.dot(M).data gseE?fastperp.dot(N).data]);')


if 0
c_eval('mvaEdJ?vec = irf.ts_vec_xyz(gseEdJ?vec.time,[gseEdJ?vec.dot(L).data gseEdJ?vec.dot(M).data gseEdJ?vec.dot(N).data]);')
c_eval('mvaEdJe?vec = irf.ts_vec_xyz(gseEdJe?vec.time,[gseEdJe?vec.dot(L).data gseEdJe?vec.dot(M).data gseEdJe?vec.dot(N).data]);')
c_eval('mvaEdJi?vec = irf.ts_vec_xyz(gseEdJi?vec.time,[gseEdJi?vec.dot(L).data gseEdJi?vec.dot(M).data gseEdJi?vec.dot(N).data]);')
c_eval('mvaRedJ?vec = irf.ts_vec_xyz(gseRedJ?vec.time,[gseRedJ?vec.dot(L).data gseRedJ?vec.dot(M).data gseRedJ?vec.dot(N).data]);')
c_eval('mvaRedJe?vec = irf.ts_vec_xyz(gseRedJe?vec.time,[gseRedJe?vec.dot(L).data gseRedJe?vec.dot(M).data gseRedJe?vec.dot(N).data]);')
c_eval('mvaRedJi?vec = irf.ts_vec_xyz(gseRedJi?vec.time,[gseRedJi?vec.dot(L).data gseRedJi?vec.dot(M).data gseRedJi?vec.dot(N).data]);')
end

mvaRotRe = irf.ts_vec_xyz(gseRotRe.time,[gseRotRe.dot(L).data gseRotRe.dot(M).data gseRotRe.dot(N).data]);
mvaGradPe = irf.ts_vec_xyz(gseGradPe.time,[gseGradPe.dot(L).data gseGradPe.dot(M).data gseGradPe.dot(N).data]);

mvaAvE = (mvaE1+mvaE2.resample(mvaE1.time)+mvaE3.resample(mvaE1.time)+mvaE4.resample(mvaE1.time))/4; 
mvaAvVe = (mvaVe1+mvaVe2.resample(mvaVe1.time)+mvaVe3.resample(mvaVe1.time)+mvaVe4.resample(mvaVe1.time))/4; 
mvaAvVeperp = (mvaVe1perp+mvaVe2perp.resample(mvaVe1perp.time)+mvaVe3perp.resample(mvaVe1perp.time)+mvaVe4perp.resample(mvaVe1perp.time))/4; 
mvaAvVi = (mvaVi1+mvaVi2.resample(mvaVi1.time)+mvaVi3.resample(mvaVi1.time)+mvaVi4.resample(mvaVi1.time))/4; 
mvaAvB = (mvaB1+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4; 
mvaAvJ = (mvaJ1+mvaJ2.resample(mvaJ1.time)+mvaJ3.resample(mvaJ1.time)+mvaJ4.resample(mvaJ1.time))/4; mvaAvJ.units = mvaJ1.units; mvaJ.name = 'avg J_fpi LMN';
mvaAvVExB = (mvaVExB1+mvaVExB2.resample(mvaVExB1.time)+mvaVExB3.resample(mvaVExB1.time)+mvaVExB4.resample(mvaVExB1.time))/4; 
mvaAvVexB = (mvaVexB1+mvaVexB2.resample(mvaVexB1.time)+mvaVexB3.resample(mvaVexB1.time)+mvaVexB4.resample(mvaVexB1.time))/4; 
mvaAvVixB = (mvaVixB1+mvaVixB2.resample(mvaVixB1.time)+mvaVixB3.resample(mvaVixB1.time)+mvaVixB4.resample(mvaVixB1.time))/4; 
mvaAvB = (mvaB1+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4; mvaAvB.name = 'avg B LMN';

mvaR0 = (mvaR1.resample(mvaR1.time)+mvaR2.resample(mvaR1.time)+mvaR3.resample(mvaR1.time)+mvaR4.resample(mvaR1.time))/4;
c_eval('mvaRR? = mvaR?-mvaR0; mvaRR? = mvaRR?.resample(irf_time(''2015-10-16T10:33:30.000Z'',''utc>epochTT'')).data;',ic)
c_eval('[mvaVe?par,mvaVe?perp,mvaVe?PA]=irf_dec_parperp(mvaB?.resample(mvaVe?),mvaVe?);',ic)
c_eval('[mvaVi?par,mvaVi?perp,mvaVi?PA]=irf_dec_parperp(mvaB?.resample(mvaVi?),mvaVi?);',ic)

% Ohm's law MVA
units = irf_units;
e = units.e; % add the minus below
mvaAvEVexB = (mvaEVexB1 + mvaEVexB2.resample(mvaEVexB1) + mvaEVexB3.resample(mvaEVexB1) + mvaEVexB4.resample(mvaEVexB1))/4;
mvaOhmGradPe = mvaGradPe/avNe.resample(mvaGradPe.time)/e*1e-9*1e-6; mvaOhmGradPe.units = 'mV/m';
mvaOhmVexB = mvaAvVexB; mvaOhmVexB.units = 'mV/m';
mvaOhmVixB = mvaAvVixB; mvaOhmVixB.units = 'mV/m';
mvaOhmJxB_a = mvaAvJ.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_a.units = 'mV/m';
mvaOhmJxB_b = mvaJcurl.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_b.units = 'mV/m';
mvaOhmJxB_c = (mvaJxB1/ne1+mvaJxB2.resample(mvaJxB1.time)/ne2+mvaJxB3.resample(mvaJxB1.time)/ne3+mvaJxB4.resample(mvaJxB1.time)/ne4)/4/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_c.units = 'mV/m'; 
mvaOhmJxB = mvaOhmJxB_c;
end

