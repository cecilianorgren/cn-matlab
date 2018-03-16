pathPaper =  '/Users/Cecilia/Dropbox (IRFU)/MMS_Nov12/';
ic=1:4;
units = irf_units;


% Convective electric fields
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

%% Make new electric field
% Reconstruct Electric field from probe voltages
baselength = 2*14.5;
multiplier = 1.5;
offsetE1z = -2.5;
offsetE2z = -3.5;
offsetE3z = 0;
offsetE4z = 0.5;
c_eval('dslE?z = irf.ts_scalar(dcv?.time,(dcv?.data(:,6)-dcv?.data(:,5))/baselength*multiplier*1e3-offsetE?z);',1:4)
c_eval('dslE?_newz = dslE?.clone(dslE?.time,[dslE?.data(:,1:2) dslE?z.data]);',1:4)

% lowpass filter E field and subtract this part
c_eval('dbcsVexB? = dbcsVe?.cross(dmpaB?.resample(dbcsVe?))*1e-3; dbcsVexB?.units = ''mV/m'';',ic)
c_eval('dbcsVixB? = dbcsVi?.cross(dmpaB?.resample(dbcsVi?))*1e-3; dbcsVixB?.units = ''mV/m'';',ic)

ffilt = 0.15;
c_eval('[dslE?par_newz,dslE?perp_newz] = irf_dec_parperp(dmpaB?,dslE?_newz); dslE?par_newz.units =''mv/m''; dslE?perp_newz.units =''mv/m'';',ic)
%c_eval('dslE?_lowf = dslE?perp_newz.filt(0,ffilt,[],3);',ic)
c_eval('dslE?_lowf = dslE?_newz.filt(0,ffilt,[],3);',ic)
c_eval('dbcsVexB?_lowf = dbcsVexB?.filt(0,ffilt,[],3);',ic)
c_eval('E?_vefit = dslE?_lowf + dbcsVexB?_lowf.resample(dslE?);',ic)
c_eval('dslE?_detrend = dslE?_newz - E?_vefit;',ic);
c_eval('dslE?_new = dslE?_detrend; dslE?_new.name = ''E new'';',ic)

c_eval('[dslE?par_new,dslE?perp_new] = irf_dec_parperp(dmpaB?,dslE?_new); dslE?par_new.name = ''E par''; dslE?perp_new.name = ''E perp'';',ic)

%[~,computername]=system('hostname');
%if strfind(computername,'ift0227887')
  load(['/Users/cno062/Data/MMS/' dirName '/defatt.mat'])
%else
%  load /Users/Cecilia/Data/MMS/2015Nov12/defatt.mat
%end

c_eval('gseE?_new = mms_dsl2gse(dslE?_new,defatt?);',ic)
c_eval('[gseE?par_new,gseE?perp_new] = irf_dec_parperp(gseB?,gseE?_new); gseE?par_new.name = ''E par''; gseE?perp_new.name = ''E perp'';',ic)

c_eval('gseEVexB?_new = gseE?_new.resample(gseVexB?.time) + gseVexB?; gseEVexB?_new.name = ''E+VexB'';',ic)
%c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

gseAvE_new = 0.25*(gseE1_new + gseE2_new.resample(gseE1_new) + gseE3_new.resample(gseE1_new) + gseE4_new.resample(gseE1_new)); gseAvE_new.name = 'av E new';
%gseAvVexB = 0.25*(gseVexB1 + gseVexB2.resample(gseVexB1) + gseVexB3.resample(gseVexB1) + gseVexB4.resample(gseVexB1)); gseVexB.name = 'av Ve new';
c_eval('gseE?_old = gseE?; gseE? = gseE?_new;',ic)

% Current from magnetic field (curlometer) 
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';

% Currents from moments
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE''; gseJe?.name = ''Je'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE''; gseJi?.name = ''Ji'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?); gseJ?.name = ''J (FPI)'';',ic);
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4;  gseAvJ.name = 'J (FPI 4sc)'; 

%% Pressure and temperature divergences
% no diagonal terms
c_eval('gsePe?_offdiag = gsePe?;',1:4)
c_eval('gsePe?_offdiag.data(:,!,!) = zeros(gsePe?_offdiag.length,1);',1:4,1:3)
c_eval('gsePe?_diag = gsePe?-gsePe?_offdiag;',1:4)
gseGradPe_offdiag = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1_offdiag,gsePe2_offdiag,gsePe3_offdiag,gsePe4_offdiag); gseGradPe_offdiag.units = 'nPa/km';  gseGradPe_offdiag.name = 'div Pe';
gseGradPe_diag = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1_diag,gsePe2_diag,gsePe3_diag,gsePe4_diag); gseGradPe_offdiag.units = 'nPa/km';  gseGradPe_offdiag.name = 'div Pe';

gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km';  gseGradPe.name = 'div Pe';
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';  gseGradPi.name = 'div Pi';
gseGradTe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';  gseGradTe.name = 'div Te';
gseGradTi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';  gseGradTi.name = 'div Ti';
gseGradNe = c_4_grad('gseR?','ne?','grad');  gseGradNe.units = 'cc/km';  gseGradNe.name = 'grad Ne';
gseGradNi = c_4_grad('gseR?','ni?','grad');  gseGradNi.units = 'cc/km';  gseGradNi.name = 'grad Ni';

avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
avTe = (gseTe1.trace/3+gseTe2.trace.resample(gseTe1.time)/3+gseTe3.trace.resample(gseTe1.time)/3+gseTe4.trace.resample(gseTe1.time)/3)/4;  avTe.name = '<Te>';
avPe = (gsePe1.trace/3+gsePe2.trace.resample(gsePe1.time)/3+gsePe3.trace.resample(gsePe1.time)/3+gsePe4.trace.resample(gsePe1.time)/3)/4;  avPe.name = '<Pe>';
avNi = (ni1+ni2.resample(ni1.time)+ni3.resample(ni1.time)+ni4.resample(ni1.time))/4; avNi.name = '<ni>';
avTi = (gseTi1.trace/3+gseTi2.trace.resample(gseTi1.time)/3+gseTi3.trace.resample(gseTi1.time)/3+gseTi4.trace.resample(gseTi1.time)/3)/4;  avTi.name = '<Ti>';
avPi = (gsePi1.trace/3+gsePi2.trace.resample(gsePi1.time)/3+gsePi3.trace.resample(gsePi1.time)/3+gsePi4.trace.resample(gsePi1.time)/3)/4;  avPi.name = '<Pi>';

epsNe = gseGradNe/avNe; epsNe.name = 'grad(ne)/ne'; epsNe.units = '1/km';
epsTe = gseGradTe/avTe; epsTe.name = 'grad(Te)/Te'; epsTe.units = '1/km';
epsPe = gseGradPe/avPe; epsPe.name = 'grad(Pe)/Pe'; epsPe.units = '1/km';
epsNi = gseGradNi/avNi; epsNi.name = 'grad(ni)/ni'; epsNi.units = '1/km';
epsTi = gseGradTi/avTi; epsTi.name = 'grad(Ti)/Ti'; epsTi.units = '1/km';
epsPi = gseGradPi/avPi; epsPi.name = 'grad(Pi)/Pi'; epsPi.units = '1/km';
epsTiNe = avTi.resample(gseGradPe)*gseGradNe/avNe; epsPi.name = 'grad(Pi)/Pi'; epsPi.units = '1/km';

%%

% Electron pressure divergence with diagonal terms set to zero
c_eval('gsePe?_offdial = gsePe?; gsePe?_offdial.data(:,1,1) = 0; gsePe?_offdial.data(:,2,2) = 0; gsePe?_offdial.data(:,3,3) = 0;',ic)
gseGradPe_offdial = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1_offdial,gsePe2_offdial,gsePe3_offdial,gsePe4_offdial); gseGradPe.units = 'nPa/km';



% ExB drift
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3; gseVExB?.units = '''';',ic) % km/s

% Magnetic field curvature 
c_eval('R? = gseR?.resample(gseB1);',1:4)
c_eval('B? = gseB?.resample(gseB1);',1:4)
[gseCurvB,avB]=c_4_grad('R?','B?','curvature'); gseCurvB.name = 'curv B'; gseCurvB.coordinateSystem = 'GSE';
curvBradius = 1/gseCurvB.abs; curvBradius.name = 'R_c';

% Magnetic shear, magnetic angle with respect to some reference time
trefshear = irf_time('2015-11-12T07:19:10.000Z','utc>epochTT');
c_eval('shearB? = irf.ts_scalar(gseB?.time,acosd(gseB?.norm.dot(gseB?.resample(trefshear).norm.data).data));',ic)


% Electron inertial term
c_eval('R? = gseR?.resample(gseB1);',1:4)
c_eval('Ve? = gseVe?.resample(gseB1);;',1:4)
[gseCurvVe,avVE]=c_4_grad('R?','Ve?','curvature'); gseCurvVe.name = 'curv Ve'; gseCurvVe.coordinateSystem = 'GSE';



% Pitchangle distribution
if 1
  try
    load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
  catch
    c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
  end
else
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
end
if 0
  %%
  c_eval('ePitch?par = ePDist?.pitchangles(dmpaB?,[0 15]);',ic)
  c_eval('ePitch?perp = ePDist?.pitchangles(dmpaB?,[75 105]);',ic)
  c_eval('ePitch?apar = ePDist?.pitchangles(dmpaB?,[165 180]);',ic)
  c_eval(['orig_size = size(ePitch?par.data);',...
          'eAparReshape = reshape(ePitch?apar.data,1, prod(orig_size));',...
          'eParReshape = reshape(ePitch?par.data,1, prod(orig_size));',...
          'eParAparReshape = eParReshape./eAparReshape;',...
          'eParApar = reshape(eParAparReshape,orig_size);',...
          'eAnis?_app = ePDist?.clone(ePitch?par.time,log(eParApar));'],ic) 
end


% T = kg?s?2?A?1
% J = kg?m2?s?2
% J/T = kg?m2?s?2/(kg?s?2?A?1) = Am2

% c_eval('Wref = 0.5*units.me*vte?perp.resample(tref).data.^2*10^6/units.eV;',ic)
% c_eval('Wref = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
% c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
%     c_eval('muref = 0.5*units.me*vte?perp.resample(tref).data.^2*10^6/units.eV/mvaB?.abs.resample(tref).data;',ic)
%     c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*muref);',ic)
    


%% Wave analysis
if 0
  %% Field aligned wavelets
  c_eval('facE? = irf_convert_fac(gseE?,gseB?,[1 0 0]);',ic)
  c_eval('facB?scm = irf_convert_fac(gseB?scm,gseB?,[1 0 0]);',ic)
  
  c_eval('wavE?fac = irf_wavelet(facE?.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavE?fac.f_units = ''Hz''; wavE?fac.f_label = ''f [Hz]''; wavE?fac.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
  c_eval('wavB?fac = irf_wavelet(facB?scm.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavB?fac.f_units = ''nT''; wavB?fac.f_label = ''f [Hz]''; wavB?fac.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
%%
if 0
  ic = 1:4;
  tic
  c_eval('wavE? = irf_wavelet(gseE?.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
  c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
    toc
%%
tintPol = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:52.00Z');
c_eval('polarization? = irf_ebsp(gseE?.tlim(tintPol),gseB?scm.tlim(tintPol),gseB?.tlim(tintPol),gseB?.tlim(tintPol),gseR?brsttime.tlim(tintPol),[10 2200],''polarization'',''fac'');',1)
end
end
%% Calculate some additional parameters
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
c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gseB?.time,rp?)*1e-3; rp?.units = ''km''; re?.name=''p gyroradius'';',ic)

%c_eval('beta?_ = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('PP? = gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3; PP?.name = ''Plasma pressure''; PB?.units =''nPa'';',ic)
c_eval('Ptot? = PB?.resample(PP?) + PP?; PP?.name = ''Plasma pressure''; PB?.units =''nPa'';',ic)
c_eval('betae? = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('betai? = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = PP?/PB?.resample(PP?);',ic)

% JxB
c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?));',ic)

% Electron pitchangles
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?);',ic);

% Non-ideal electric field, E+VexB
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

% Perpendicular and parallel velocity and current components
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)

% Parallel and perpendicular electric fields
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
try
c_eval('[gseE?hmfepar,gseE?hmfeperp] = irf_dec_parperp(gseB?,gseE?hmfe); gseE?hmfepar.name = ''E par hmfe''; gseE?hmfeperp.name = ''E perp hmfe'';',ic)
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

c_eval('gseVexBE? = gseE?.resample(gseVe?)+gseVexB?.data;',ic);

c_eval('deltae? = gseVexB?/(gseVe?perp.abs*gseB?.resample(gseVe?).abs*1e-9);',ic);
%c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

gseR0 = (gseR1.resample(gseR1.time)+gseR2.resample(gseR1.time)+gseR3.resample(gseR1.time)+gseR4.resample(gseR1.time))/4;
c_eval('gseRR? = gseR?-gseR0; gseRR? = gseRR?.resample(irf_time(''2015-11-12T07:19:21.000Z'',''utc>epochTT'')).data;',ic)


% Plasma beta and magnetic pressure
c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

c_eval('wavVe?par = irf_wavelet(gseVe?par.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?par.f_units = ''Hz''; wavVe?par.f_label = ''f [Hz]''; wavVe?par.p_label = {''log_{10} v_{e,||}^2'',''(km/s)^2/Hz''};',ic)
c_eval('wavVe?perp = irf_wavelet(gseVe?perp.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?perp.f_units = ''Hz''; wavVe?perp.f_label = ''f [Hz]''; wavVe?perp.p_label = {''log_{10} v_{e,\perp}^2'',''(km/s)^2/Hz''};',ic)


%% Average properties

gseAvE = (gseE1+gseE2.resample(gseE1.time)+gseE3.resample(gseE1.time)+gseE4.resample(gseE1.time))/4; gseAvE.name = 'avg E (gse)';
gseAvE_old = (gseE1_old+gseE2_old.resample(gseE1.time)+gseE3_old.resample(gseE1.time)+gseE4_old.resample(gseE1.time))/4; gseAvE_old.name = 'avg E (gse)';

avPB = (PB1 + PB2.resample(PB1) + PB3.resample(PB1) + PB4.resample(PB1))/4;

gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; 
gseAvVeperp = (gseVe1perp+gseVe2perp.resample(gseVe1perp.time)+gseVe3perp.resample(gseVe1perp.time)+gseVe4perp.resample(gseVe1perp.time))/4; 
gseAvVi = (gseVi1+gseVi2.resample(gseVi1.time)+gseVi3.resample(gseVi1.time)+gseVi4.resample(gseVi1.time))/4; 
avTi = (gseTi1.trace/3+gseTi2.trace.resample(gseTi1.time)/3+gseTi3.trace.resample(gseTi1.time)/3+gseTi4.trace.resample(gseTi1.time)/3)/4; 
avPi = (gsePi1.trace/3+gsePi2.trace.resample(gsePi1.time)/3+gsePi3.trace.resample(gsePi1.time)/3+gsePi4.trace.resample(gsePi1.time)/3)/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; 
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.units = gseJ1.units;
gseAvVExB = (gseVExB1+gseVExB2.resample(gseVExB1.time)+gseVExB3.resample(gseVExB1.time)+gseVExB4.resample(gseVExB1.time))/4; 
gseAvEVexB = (gseEVexB1+gseEVexB2.resample(gseEVexB1.time)+gseEVexB3.resample(gseEVexB1.time)+gseEVexB4.resample(gseEVexB1.time))/4; 
gseAvVexB = (gseVexB1+gseVexB2.resample(gseVexB1.time)+gseVexB3.resample(gseVexB1.time)+gseVexB4.resample(gseVexB1.time))/4; 
gseAvVixB = (gseVixB1+gseVixB2.resample(gseVixB1.time)+gseVixB3.resample(gseVixB1.time)+gseVixB4.resample(gseVixB1.time))/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4;

[gseGradPepar,gseGradPeperp] = irf_dec_parperp(gseAvB,gseGradPe); gseGradPepar.name = 'grad Pe par'; gseGradPeperp.name = 'grad Pe perp';
% Inertial term
divVe = c_4_grad('gseR?','gseVe?','div'); divVe.name = 'div Ve';

% Ohm's law terms
%gseAvEVexB = (gseEVexB1 + gseEVexB2.resample(gseEVexB1) + gseEVexB3.resample(gseEVexB1) + gseEVexB4.resample(gseEVexB1))/4;
gseOhmGradPe = gseGradPe/avNe.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseOhmGradPe.units = 'mV/m';
gseOhmVexB = gseAvVexB; gseOhmVexB.units = 'mV/m';
gseOhmVixB = gseAvVixB; gseOhmVixB.units = 'mV/m';
gseOhmJxB_a = gseAvJ.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_a.units = 'mV/m';
gseOhmJxB_b = gseJcurl.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_b.units = 'mV/m';
gseOhmJxB_c = (gseJxB1/ne1+gseJxB2.resample(gseJxB1.time)/ne2+gseJxB3.resample(gseJxB1.time)/ne3+gseJxB4.resample(gseJxB1.time)/ne4)/4/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_c.units = 'mV/m'; 
gseOhmJxB = gseOhmJxB_c;



% Other parameters 
[gseGradPepar,gseGradPeperp] = irf_dec_parperp(gseAvB,gseGradPe); gseGradPepar.name = 'div Pe par'; gseGradPeperp.name = 'div Pe par';
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

% Rotation of E
[gseRotE,E]=c_4_grad('gseR?','gseE?','curl');

% Magnetic moment
c_eval('mag_mom? = 0.5*units.me*vte?perp.^2*10^6/(gseB?.abs*1e-9)*1e9;  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)
%% MVA: Rotate into lmn coordinates
ic=1:4;
% Set up coordinate system
% Default
%[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
%L = v(1,:); M = v(2,:); N = v(3,:);
%coordLabels = {'L','M','N'};
%lmn = [N;-M;L];


%% compare mva's, large scale
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(irf.tint(''2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z'')));',1:4)
%v2 = [v2(1,:); -v2(2,:); -v2(3,:)];
clear mva_angle1 mva_angle2 mva_angle3
mva_angle1 = []; mva_angle2 = []; mva_angle3 = [];
for ic1 = 1:3
  for ic2 = ic1:4
    for icomp = 1:3
      c_eval('mva_angle_temp = real(acosd(dot(v?(icomp,:),v!(icomp,:))));',ic1,ic2)
      c_eval('mva_angle? = [mva_angle? mva_angle_temp];',icomp)
    end
  end
end
c_eval('lratio?large(1) = l?(1)/l?(2); lratio?large(2) = l?(2)/l?(3);',1:4)
lratio_mean_large = mean([lratio1large; lratio2large; lratio3large; lratio4large]);
lratio_std_large = std([lratio1large; lratio2large; lratio3large; lratio4large]);


c_eval('mva_all(:,:,?) = v?;',1:4)
mva_mean = mean(mva_all,3);
mva_angle_mean = mean([mva_angle1' mva_angle2' mva_angle3'],1);
mva_angle_std = std([mva_angle1' mva_angle2' mva_angle3'],1);

mva_mean_large = mva_mean;
mva_angle_mean_large = mva_angle_mean;
mva_angle_std_large = mva_angle_std;

%% compare mva's, befurcated current sheet, enforced normal component
tint_bcs_utc = '2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z';
tint_bcs = irf.tint(tint_bcs_utc);
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(tint_bcs),''td'');',1:4)
c_eval('v? = [v?(1,:); -v?(2,:); -v?(3,:)];',[1 2 3 4])
clear mva_angle1 mva_angle2 mva_angle3
mva_angle1 = []; mva_angle2 = []; mva_angle3 = [];
for ic1 = 1:3
  for ic2 = ic1:4
    
    for icomp = 1:3
      c_eval('mva_angle_temp = real(acosd(dot(v?(icomp,:),v!(icomp,:))));',ic1,ic2)
      c_eval('mva_angle? = [mva_angle? mva_angle_temp];',icomp)
    end
  end
end
c_eval('lratio?(1) = l?(1)/l?(2); lratio?(2) = l?(2)/l?(3);',1:4)
lratio_mean_bcs_td = mean([lratio1; lratio2; lratio3; lratio4]);
lratio_std_bcs_td = std([lratio1; lratio2; lratio3; lratio4]);

c_eval('mva_all_td(:,:,?) = v?;',1:4)
mva_mean_td = mean(mva_all_td,3);
mva_angle_mean_td = mean([mva_angle1' mva_angle2' mva_angle3'],1);
mva_angle_std_td = std([mva_angle1' mva_angle2' mva_angle3'],1);

c_eval('mva_angle_large_bifurcated_td(?) = real(acosd(dot(mva_mean_td(?,:),mva_mean_large(?,:))));',1:3)

%% compare mva's, befurcated current sheet
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(tint_bcs));',1:4)
v2 = [v2(1,:); -v2(2,:); -v2(3,:)];
clear mva_angle1 mva_angle2 mva_angle3
mva_angle1 = []; mva_angle2 = []; mva_angle3 = [];
for ic1 = 1:3
  for ic2 = ic1:4
    
    for icomp = 1:3
      c_eval('mva_angle_temp = real(acosd(dot(v?(icomp,:),v!(icomp,:))));',ic1,ic2)
      c_eval('mva_angle? = [mva_angle? mva_angle_temp];',icomp)
    end
  end
end
c_eval('lratio?(1) = l?(1)/l?(2); lratio?(2) = l?(2)/l?(3);',1:4)
lratio_mean_bcs = mean([lratio1; lratio2; lratio3; lratio4]);
lratio_std_bcs = std([lratio1; lratio2; lratio3; lratio4]);

c_eval('mva_all(:,:,?) = v?;',1:4)
mva_mean = mean(mva_all,3);
mva_angle_mean = mean([mva_angle1' mva_angle2' mva_angle3'],1);
mva_angle_std = std([mva_angle1' mva_angle2' mva_angle3'],1);

c_eval('mva_angle_large_bifurcated(?) = real(acosd(dot(mva_mean(?,:),mva_mean_large(?,:))));',1:3)

%% compare mva's: large and small cs after the bifurcated one
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(irf.tint(''2015-11-12T07:19:30.200Z/2015-11-12T07:19:31.200Z'')));',1:4)
c_eval('v? = [v?(1,:); -v?(2,:); -v?(3,:)];',1:4)
clear mva_angle1 mva_angle2 mva_angle3
mva_angle1 = []; mva_angle2 = []; mva_angle3 = [];
for ic1 = 1:3
  for ic2 = ic1:4
    for icomp = 1:3
      c_eval('mva_angle_temp = real(acosd(dot(v?(icomp,:),v!(icomp,:))));',ic1,ic2)
      c_eval('mva_angle? = [mva_angle? mva_angle_temp];',icomp)
    end
  end
end
c_eval('mva_all(:,:,?) = v?;',1:4)
mva_mean_2 = mean(mva_all,3);
mva_angle_mean_2 = mean([mva_angle1' mva_angle2' mva_angle3'],1);
mva_angle_std_2 = std([mva_angle1' mva_angle2' mva_angle3'],1);


c_eval('mva_angle_bcs_other_small(?) = real(acosd(dot(mva_mean(?,:),mva_mean_2(?,:))));',1:3)

if 0
  %%
  figure(17)
  plot_quivers(mva_mean,mva_mean*0,'b',{'L','M','N'}); hold on
  plot_quivers(mva_mean_2,mva_mean*0,'r',{'L','M','N'});
  plot_quivers(mva_mean_large,mva_mean*0,'k',{'L','M','N'});
  set(gca,'ylim',[-1 1],'xlim',[-1 1],'zlim',[-1 1])
  
end

%% Choose LMN coordinate system
coordSystem = 23; % 16:td, 15:unconstr.
switch coordSystem % Choose another
  case 1 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
    L = v(1,:); M = -v(2,:); N = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 12 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 15 % N: minimum variance of B, mean of all 4 spacecraft
    %[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z')));
    L = mva_mean(1,:); M = mva_mean(2,:); N = mva_mean(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];
  case 16 % N: minimum variance of B, mean of all 4 spacecraft
    %[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z')),'<Bn>=0');
    L = mva_mean_td(1,:); M = mva_mean_td(2,:); N = mva_mean_td(3,:);
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
    [out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
    L = -v(2,:); M = -v(1,:); N = v(3,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
    [out,l,v] = irf_minvar(gseJcurl.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];   
  case 23 % N: minimum variance of J, L: maximum variance of B
    %%
    tint_minvar = irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z');
    L_B = mva_mean(1,:); M_B = mva_mean(2,:); N_B = mva_mean(3,:); 
    %mva_mean
    if 0
    [outB,lB,vB] = irf_minvar(gseAvB.tlim(tint_minvar));
    %vB
    L_B = vB(1,:); M_B = -vB(2,:); N_B = -vB(3,:);
    vB = [L_B;M_B;N_B];
    end        
    
    [outJ,lJ,vJ] = irf_minvar(gseJcurl.tlim(tint_minvar));    
    L_J = vJ(1,:); M_J = vJ(2,:); N_J = vJ(3,:);
    vJ = [L_J;M_J;N_J];
    if 1
      L = L_B;
      N = cross(L,cross(N_J,L));
      M = cross(N,L);
    elseif 0      
      N = N_J;
      L = cross(N,cross(L_B,N));
      M = cross(N,L);
    else
      N = N_B;
      L = L_B;
      M = M_B;
    end
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

% DeHoffman-Teller frame
[vht,eht,dvht,p,cc]=irf_vht(-gseVixB1.tlim(tint_bcs+[-2 3]),gseB1.tlim(tint_bcs+[-2 3]));
mvaVht = vht*lmn';
mvaEht = eht*lmn';

disp(sprintf('L = [%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L,M,N))
% Rotate data
c_eval('mvaE?_new = gseE?_new*lmn'';',ic)
c_eval('mvaE?par_new = gseE?par_new; mvaE?perp_new = gseE?perp_new*lmn'';',ic)
mvaAvE_new = gseAvE_new*lmn'';
c_eval('mvaEVexB?_new = gseEVexB?_new*lmn'';',ic)

c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';')
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';')
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';')
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';')
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';')
c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';')
c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';')
c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';')
mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)

c_eval('mvaVexB? = irf.ts_vec_xyz(gseVexB?.time,[gseVexB?.dot(L).data gseVexB?.dot(M).data gseVexB?.dot(N).data]); mvaVexB?.units = ''mV/m'';')
c_eval('mvaVixB? =  irf.ts_vec_xyz(gseVixB?.time,[gseVixB?.dot(L).data gseVixB?.dot(M).data gseVixB?.dot(N).data]); mvaVixB?.units = ''mV/m'';')
c_eval('mvaEVexB? =  irf.ts_vec_xyz(gseEVexB?.time,[gseEVexB?.dot(L).data gseEVexB?.dot(M).data gseEVexB?.dot(N).data]); mvaEVexB?.units = ''mV/m'';')
%mvaVDe =  irf.ts_vec_xyz(vDe.time,[vDe.dot(L).data vDe.dot(M).data vDe.dot(N).data]); mvaVDe.units = '';
mvaAvJ =  gseAvJ*lmn'; mvaAvJ.units = 'nA/m^2';
c_eval('mvaJxB? = mvaJ?.cross(mvaB?.resample(mvaJ?.time));')
c_eval('mvaVExB? =  gseVExB?*lmn'';')
c_eval('mvaVe?par = gseVe?par;')
c_eval('mvaVe?perp = gseVe?perp*lmn''; mvaVe?perp.name = ''Ve perp lmn'';')
c_eval('mvaJ?par = gseJ?par; mvaJ?par.units = ''nA/m^2'';')
c_eval('mvaJ?perp = gseJ?perp*lmn''; mvaJ?perp.name = ''J perp lmn''; mvaJ?perp.units = ''nA/m^2'';')
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

mvaRotE = gseRotE*lmn';

mvaRotRe = irf.ts_vec_xyz(gseRotRe.time,[gseRotRe.dot(L).data gseRotRe.dot(M).data gseRotRe.dot(N).data]);
mvaGradPi = irf.ts_vec_xyz(gseGradPi.time,[gseGradPi.dot(L).data gseGradPi.dot(M).data gseGradPi.dot(N).data]);
mvaGradPe = irf.ts_vec_xyz(gseGradPe.time,[gseGradPe.dot(L).data gseGradPe.dot(M).data gseGradPe.dot(N).data]);
mvaGradPe_diag = irf.ts_vec_xyz(gseGradPe_diag.time,[gseGradPe_diag.dot(L).data gseGradPe_diag.dot(M).data gseGradPe_diag.dot(N).data]);
mvaGradPe_offdiag = irf.ts_vec_xyz(gseGradPe_offdiag.time,[gseGradPe_offdiag.dot(L).data gseGradPe_offdiag.dot(M).data gseGradPe_offdiag.dot(N).data]);

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
c_eval('mvaRR? = mvaR?-mvaR0; mvaRR? = mvaRR?.resample(irf_time(''2015-11-12T07:19:21.000Z'',''utc>epochTT'')).data;',ic)
c_eval('[mvaVe?par,mvaVe?perp,mvaVe?PA]=irf_dec_parperp(mvaB?.resample(mvaVe?),mvaVe?);',ic)
c_eval('[mvaVi?par,mvaVi?perp,mvaVi?PA]=irf_dec_parperp(mvaB?.resample(mvaVi?),mvaVi?);',ic)

% Ohm's law MVA
units = irf_units;
e = units.e; % add the minus below
mvaAvEVexB = (mvaEVexB1 + mvaEVexB2.resample(mvaEVexB1) + mvaEVexB3.resample(mvaEVexB1) + mvaEVexB4.resample(mvaEVexB1))/4;
mvaOhmGradPe = mvaGradPe/avNe.resample(mvaGradPe.time)/e*1e-9*1e-6; mvaOhmGradPe.units = 'mV/m';
mvaOhmGradPi = mvaGradPi/avNe.resample(mvaGradPi.time)/e*1e-9*1e-6; mvaOhmGradPi.units = 'mV/m';
mvaOhmGradPe_diag = mvaGradPe_diag/avNe.resample(mvaGradPe.time)/e*1e-9*1e-6; mvaOhmGradPe_diag.units = 'mV/m';
mvaOhmGradPe_offdiag = mvaGradPe_offdiag/avNe.resample(mvaGradPe.time)/e*1e-9*1e-6; mvaOhmGradPe_offdiag.units = 'mV/m';
mvaOhmVexB = mvaAvVexB; mvaOhmVexB.units = 'mV/m';
mvaOhmVixB = mvaAvVixB; mvaOhmVixB.units = 'mV/m';
c_eval('mvaOhmJxB? = mvaJ?.resample(mvaB?.time).cross(mvaB?)/ne?/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB?.units = ''mV/m'';',1:4)
mvaOhmJxB_a = mvaAvJ.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_a.units = 'mV/m';
mvaOhmJxB_b = mvaJcurl.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_b.units = 'mV/m';
mvaOhmJxB_c = (mvaJxB1/ne1+mvaJxB2.resample(mvaJxB1.time)/ne2+mvaJxB3.resample(mvaJxB1.time)/ne3+mvaJxB4.resample(mvaJxB1.time)/ne4)/4/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_c.units = 'mV/m'; 
mvaOhmJxB = mvaOhmJxB_c;

% Curvature
c_eval('R? = mvaR?.resample(mvaB1);',1:4)
c_eval('B? = mvaB?.resample(mvaB1);',1:4)
[mvaCurvB,BB]=c_4_grad('R?','B?','curvature'); mvaCurvB.name = 'curv B'; mvaCurvB.units = '1/km';
curvBradius = 1/mvaCurvB.abs; curvBradius.name = 'R_c';


figure(12)
c_eval('h=irf_plot({mvaB?,mvaJ?,mvaVe?});',1)
%irf_zoom(h,'x',tintZoom); irf_zoom(h,'y');
h(1).Title.String = sprintf('L = [%.2f %.2f %.2f]; M = [%.2f %.2f %.2f]; N = [%.2f %.2f %.2f];',L,M,N);
%fprintf('L = [%.2f %.2f %.2f]; M = [%.2f %.2f %.2f]; N = [%.2f %.2f %.2f];\n',lmn')
%fprintf('L = [%.2f %.2f %.2f]; M = [%.2f %.2f %.2f]; N = [%.2f %.2f %.2f];\n',L,M,N)

%%

CS_normal_velocity = 70;
c_eval('dtP = mvaPe?.time(2)-mvaPe?.time(1); LPNN?time = mvaPe?.time(1:end-1)+0.5*dtP;',ic)
c_eval('sumP? = mvaPe?.zz.data(1:end-1)+mvaPe?.zz.data(2:end);',ic)
c_eval('LPNN? = irf.ts_scalar(LPNN?time,0.5*sumP?./diff(mvaPe?.zz.data(:,1))*CS_normal_velocity*dtP);',ic)
0.5*(0.055+0.03)/(0.055-0.03)*70*dtP;
%%



c_eval('R? = gseR?.resample(ne?);',1:4)
c_eval('n? = ne?.resample(ne?);',1:4)
[gseGradNe,NN]=c_4_grad('R?','n?','grad'); gseGradNe.name = 'grad ne'; gseGradNe.units = '1/km';
gseEfromdivPi = units.kB*avTi.resample(gseGradNe)*gseGradNe/units.e*1e3; gseEfromdivPi.units = 'mV/m'; gseEfromdivPi.name = 'E (divPi)';


c_eval('R? = gseR?.resample(ne?);',1:4)
c_eval('Pi? = gsePi?.resample(ne?).trace/3;',1:4)
[gseGradPi,NN]=c_4_grad('R?','Pi?','grad'); gseGradPi.name = 'grad Pi'; gseGradPi.units = '1/km';
gseEfromdivPi_ = avPi.resample(gseGradPi)*gseGradPi/units.e*1e3; gseEfromdivPi.units = 'mV/m'; gseEfromdivPi.name = 'E (divPi)';


% c_eval('R? = mvaR?.resample(facPe1);',1:4)
% c_eval('facPdiff? = facPe?.xx-0.5*(facPe?.yy+facPe?.zz); facPdiff? = facPdiff?.resample(facPe1);',1:4)
% [mvaCurvPdiff,BB]=c_4_grad('R?','facPdiff?','curvature'); mvaCurvPdiff.name = 'curv ppar-pperp'; mvaCurvPdiff.units = '1/km';
% curvPdiffradius = 1/mvaCurvPdiff.abs; curvPdiffradius.name = 'R_c (ppar-pperp)';

c_eval('R? = mvaR?.resample(mvaVe1);',1:4)
c_eval('Ve? = mvaVe?.resample(mvaVe1);',1:4)
[mvaCurvVe,VVE]=c_4_grad('R?','Ve?','curvature'); mvaCurvVe.name = 'curv Ve';



