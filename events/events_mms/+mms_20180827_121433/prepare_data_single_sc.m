  ic = 1:3;
units = irf_units;

%% Current
% Currents from moments, use ne also for Ji 
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

%% Perpendicular and parallel decomposition
% Celocity and current
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)
c_eval('[gsmVe?par,gsmVe?perp] = irf_dec_parperp(gsmB?,gsmVe?); gsmVe?par.name = ''Ve par''; gsmVe?perp.name = ''Ve perp'';',ic)
c_eval('[gsmVi?par,gsmVi?perp] = irf_dec_parperp(gsmB?,gsmVi?); gsmVi?par.name = ''Vi par''; gsmVi?perp.name = ''Vi perp'';',ic)


% Electric fields
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('[gsmE?par,gsmE?perp] = irf_dec_parperp(gsmB?,gsmE?); gsmE?par.name = ''E par''; gsmE?perp.name = ''E perp'';',ic)
% Wave magnetic field
%gseB2scm = gseB2scm{2};
c_eval('[gseB?scmpar,gseB?scmperp] = irf_dec_parperp(gseB?,gseB?scm); gseB?scmpar.name = ''B par scm''; gseB?scmperp.name = ''B perp scm'';',ic)
try
c_eval('[gseE?hmfepar,gseE?hmfeperp] = irf_dec_parperp(gseB?,gseE?hmfe); gseE?hmfepar.name = ''E par hmfe''; gseE?hmfeperp.name = ''E perp hmfe'';',ic)
end

%% Cross products
% ExB drift
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gsmVExB? = c_coord_trans(''GSE'',''GSM'',gseVExB?);',ic)

% Convective electric fields
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)
c_eval('gsmVexB? = c_coord_trans(''GSE'',''GSM'',gseVexB?);',ic)
c_eval('gsmVixB? = c_coord_trans(''GSE'',''GSM'',gseVixB?);',ic)

% Non-ideal electric field, E+VexB
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)
c_eval('gseEVixB? = gseE?.resample(gseVixB?.time)+gseVixB?; gseEVixB?.name = ''E+VixB'';',ic)
c_eval('gsmEVexB? = c_coord_trans(''GSE'',''GSM'',gseEVexB?);',ic)
c_eval('gsmEVixB? = c_coord_trans(''GSE'',''GSM'',gseEVixB?);',ic)

% JxB
c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?));',ic)

% Magnetic field curvature 
% if all(ic==[1:4])
% c_eval('R? = gseR?.resample(gseB1);',1:4)
% c_eval('B? = gseB?.resample(gseB1);',1:4)
% [gseCurvB,avB]=c_4_grad('R?','B?','curvature'); gseCurvB.name = 'curv B'; gseCurvB.coordinateSystem = 'GSE';
% curvBradius = 1/gseCurvB.abs; curvBradius.name = 'R_c';
% end
%% Pitchangle distributions
if 0
  load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
  %load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
elseif 0
  %%
  ic = 1;
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
  c_eval('ePitch?par = ePDist?.pitchangles(dmpaB?,[0 15]);',ic)
  c_eval('ePitch?perp = ePDist?.pitchangles(dmpaB?,[75 105]);',ic)
  c_eval('ePitch?apar = ePDist?.pitchangles(dmpaB?,[165 180]);',ic)
  ic = 1:4;
end

%% Calculate some additional parameters, irf_plasma_calc
% Speeds
c_eval('matB? = gseB?.abs.data;',ic)
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
c_eval('flh? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Flh''); flh? = irf.ts_scalar(gseB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fce''); fce? = irf.ts_scalar(gseB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fcp''); fcp? = irf.ts_scalar(gseB?.time,fcp?);',ic)
c_eval('fpe? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpe''); fpe? = irf.ts_scalar(gseB?.time,fpe?);',ic)
c_eval('fpp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpp''); fpp? = irf.ts_scalar(gseB?.time,fpp?);',ic)

% Length scales
c_eval('Lp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Li''); Lp? = irf.ts_scalar(gseB?.time,Lp?)*1e-3; Lp?.units = ''km''; Lp?.name=''p inertial length'';',ic)
c_eval('Le? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Le''); Le? = irf.ts_scalar(gseB?.time,Le?)*1e-3; Le?.units = ''km''; Le?.name=''e inertial length'';',ic)
c_eval('Ld? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Ld''); Ld? = irf.ts_scalar(gseB?.time,Ld?)*1e-3; Ld?.units = ''km''; Ld?.name=''Debye length'';',ic)
c_eval('re? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Roe''); re? = irf.ts_scalar(gseB?.time,re?)*1e-3; re?.units = ''km''; re?.name=''e gyroradius'';',ic)
c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gseB?.time,rp?)*1e-3; rp?.units = ''km''; rp?.name=''p gyroradius'';',ic)

%c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

% Magnetic moment
c_eval('mag_mom? = 0.5*units.me*vte?perp.^2*10^6/(gseB?.abs*1e-9)*1e9;  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)

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

% Plasma beta and magnetic pressure
%c_eval('beta?_ = (re?/Le?).^2;',ic) % this is beta_e
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

c_eval('wavVe?par = irf_wavelet(gseVe?par.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?par.f_units = ''Hz''; wavVe?par.f_label = ''f [Hz]''; wavVe?par.p_label = {''log_{10} v_{e,||}^2'',''(km/s)^2/Hz''};',ic)
c_eval('wavVe?perp = irf_wavelet(gseVe?perp.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?perp.f_units = ''Hz''; wavVe?perp.f_label = ''f [Hz]''; wavVe?perp.p_label = {''log_{10} v_{e,\perp}^2'',''(km/s)^2/Hz''};',ic)

%% LMN
tint_zoom = irf.tint('2018-08-27T12:15:35Z/2018-08-27T12:15:50Z');

L = [0.91, -0.41, 0.10];
M = [0.41, 0.91, -0.02];
N = [-0.08, 0.06, 0.99];
lmn = [L;M;N];
% try tiny rotation in MN plane
angle = 14;
Rmat = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
MNnew = Rmat*[M;N];
lmn_new = [L;MNnew];
N = MNnew(2,:);
M = MNnew(1,:);
lmn = lmn_new;

disp(sprintf('L = [%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L,M,N))
% Rotate data
% c_eval('mvaR? = gsmR?*lmn''; mvaR?.name = ''R LMN'';',ic)
c_eval('mvaB? = gsmB?*lmn''; mvaB?.name = ''B LMN'';',ic)
%c_eval('mvaB?scm = gsmB?scm*lmn''; mvaB?scm.name = ''B LMN'';')
c_eval('mvaE? = gsmE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gsmVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gsmVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
% c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';')
% c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';')
% c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';')
% mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
% c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
% c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
% c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
% c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)

c_eval('mvaVexB? = irf.ts_vec_xyz(gsmVexB?.time,[gsmVexB?.dot(L).data gsmVexB?.dot(M).data gsmVexB?.dot(N).data]); mvaVexB?.units = ''mV/m'';',ic)
c_eval('mvaVixB? =  irf.ts_vec_xyz(gsmVixB?.time,[gsmVixB?.dot(L).data gsmVixB?.dot(M).data gsmVixB?.dot(N).data]); mvaVixB?.units = ''mV/m'';',ic)
c_eval('mvaEVexB? =  irf.ts_vec_xyz(gsmEVexB?.time,[gsmEVexB?.dot(L).data gsmEVexB?.dot(M).data gsmEVexB?.dot(N).data]); mvaEVexB?.units = ''mV/m'';',ic)
%mvaVDe =  irf.ts_vec_xyz(vDe.time,[vDe.dot(L).data vDe.dot(M).data vDe.dot(N).data]); mvaVDe.units = '';
%mvaAvJ =  gseAvJ*lmn'; mvaAvJ.units = 'nA/m^2';
%c_eval('mvaJxB? = mvaJ?.cross(mvaB?.resample(mvaJ?.time));',ic)
c_eval('mvaVExB? =  gsmVExB?*lmn'';',ic)
c_eval('mvaVe?par = gsmVe?par;',ic)
c_eval('mvaVe?perp = gsmVe?perp*lmn''; mvaVe?perp.name = ''Ve perp lmn'';',ic)
%c_eval('mvaJ?par = gseJ?par; mvaJ?par.units = ''nA/m^2'';')
%c_eval('mvaJ?perp = gseJ?perp*lmn''; mvaJ?perp.name = ''J perp lmn''; mvaJ?perp.units = ''nA/m^2'';',ic)
c_eval('mvaE?par = gseE?par;',ic)
c_eval('mvaE?perp = gseE?perp*lmn''; mvaE?perp.name = ''E perp lmn'';',ic)
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

%%
disp('Done preparing data.')