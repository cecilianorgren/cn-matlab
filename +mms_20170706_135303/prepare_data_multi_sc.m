ic = 1:4;
units = irf_units;

%% Tetrahedron centered positions 9
if all(ic==[1:4])
gseR0 = (gseR1.resample(gseR1.time)+gseR2.resample(gseR1.time)+gseR3.resample(gseR1.time)+gseR4.resample(gseR1.time))/4;
c_eval('gseRR? = gseR?-gseR0; gseRR? = gseRR?.resample(irf_time(''2017-07-11T22:33:28.00Z'',''utc>epochTT'')).data;',ic)
end

%% Current, curlometer
if all(ic==[1:4])
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
end

%% Pressure and temperature divergences
if all(ic==[1:4])
gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km'; gseGradPe.name = 'div Pe';
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';
gseGradTe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';
gseGradTi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';
gseGradNe = c_4_grad('gseR?','ne?','grad');
gseGradNi = c_4_grad('gseR?','ni?','grad');
end

%% Magnetic field curvature 
if all(ic==[1:4])
c_eval('R? = gseR?.resample(gseB1);',1:4)
c_eval('B? = gseB?.resample(gseB1);',1:4)
[gseCurvB,avB]=c_4_grad('R?','B?','curvature'); gseCurvB.name = 'curv B'; gseCurvB.coordinateSystem = 'GSE';
curvBradius = 1/gseCurvB.abs; curvBradius.name = 'R_c';
end

%% Divergence of E, charge pertubation
if 1
  %%
  c_eval('R? = gseR?.resample(gseE1)*1e3;',1:4)
  c_eval('E? = gseE?.resample(gseE1)*1e-3;',1:4)
  [divE,avE]=c_4_grad('R?','E?','div'); divE.name = 'div E';
  dn_divE = divE*units.eps0/units.e*1e-6; 
  dn_divE.name = 'dn from div E';
  dn_divE.units = 'cm^-3';
  %% using only Epar
  c_eval('R? = gseR?.resample(gseE1)*1e3;',1:4)
  c_eval('E?parvec = gseE?.dot(gseB?.resample(gseE?).norm)*gseB?.resample(gseE?).norm;',1:4)
  c_eval('E?_paronly =E?parvec.resample(gseE1)*1e-3;',1:4)
  [divE_paronly,avE]=c_4_grad('R?','E?_paronly','div'); divE_paronly.name = 'div E';
  dn_divE_par = divE_paronly*units.eps0/units.e*1e-6; 
  dn_divE_par.name = 'dn from div E';
  dn_divE_par.units = 'cm^-3';
end
%% Inertial term, div(ve)
if all(ic==[1:4])
  divVe = c_4_grad('gseR?','gseVe?','div'); divVe.name = 'div Ve';
end
%% Average properties
if all(ic==[1:4])
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 
avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
gseAvE = (gseE1+gseE2.resample(gseE1.time)+gseE3.resample(gseE1.time)+gseE4.resample(gseE1.time))  /4; 
gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; 
gseAvVeperp = (gseVe1perp+gseVe2perp.resample(gseVe1perp.time)+gseVe3perp.resample(gseVe1perp.time)+gseVe4perp.resample(gseVe1perp.time))/4; 
gseAvVi = (gseVi1+gseVi2.resample(gseVi1.time)+gseVi3.resample(gseVi1.time)+gseVi4.resample(gseVi1.time))/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; 
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.units = gseJ1.units;
gseAvVExB = (gseVExB1+gseVExB2.resample(gseVExB1.time)+gseVExB3.resample(gseVExB1.time)+gseVExB4.resample(gseVExB1.time))/4; 
gseAvVexB = (gseVexB1+gseVexB2.resample(gseVexB1.time)+gseVexB3.resample(gseVexB1.time)+gseVexB4.resample(gseVexB1.time))/4; 
gseAvVixB = (gseVixB1+gseVixB2.resample(gseVixB1.time)+gseVixB3.resample(gseVixB1.time)+gseVixB4.resample(gseVixB1.time))/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4;

[gseGradPepar,gseGradPeperp] = irf_dec_parperp(gseAvB,gseGradPe); gseGradPepar.name = 'grad Pe par'; gseGradPeperp.name = 'grad Pe perp';
end
%% Ohm's law terms
if all(ic==[1:4])
%gseAvEVexB = (gseEVexB1 + gseEVexB2.resample(gseEVexB1) + gseEVexB3.resample(gseEVexB1) + gseEVexB4.resample(gseEVexB1))/4;
gseOhmGradPe = gseGradPe/avNe.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseOhmGradPe.units = 'mV/m';
gseOhmVexB = gseAvVexB; gseOhmVexB.units = 'mV/m';
gseOhmVixB = gseAvVixB; gseOhmVixB.units = 'mV/m';
gseOhmJxB_a = gseAvJ.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_a.units = 'mV/m';
gseOhmJxB_b = gseJcurl.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_b.units = 'mV/m';
gseOhmJxB_c = (gseJxB1/ne1+gseJxB2.resample(gseJxB1.time)/ne2+gseJxB3.resample(gseJxB1.time)/ne3+gseJxB4.resample(gseJxB1.time)/ne4)/4/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_c.units = 'mV/m'; 
gseOhmJxB = gseOhmJxB_c;
end

%% Other parameters 
if all(ic==[1:4])
[gseGradPepar,gseGradPeperp] = irf_dec_parperp(gseAvB,gseGradPe); gseGradPepar.name = 'div Pe par'; gseGradPeperp.name = 'div Pe par';
avPe = (gsePe1.trace+gsePe2.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1))/3/4;
LgradP = avPe/gseGradPe.abs; LgradP.name = 'L_{P}';
end

% Proxy for magnetic field slippage (Scudder2105)
if all(ic==[1:4])
avFce = (fce1+fce2.resample(fce1)+fce3.resample(fce1)+fce4.resample(fce1))/4;
avRhoe = (re1+re2.resample(re1)+re3.resample(re1)+re4.resample(re1))/4;
Y = avFce*(avRhoe/LgradP).^2; Y.units = 's^-1';
Ynodim = (avRhoe/LgradP); Y.units = 's^-1';

facAvTe = (facTe1+facTe2.resample(facTe1)+facTe3.resample(facTe1)+facTe4.resample(facTe1))/4;

% Curl of E+VexB (Scudder2105)
gseRotRe = mms_2015Oct16.rotRe(gseR1,gseR2,gseR3,gseR4,gseEVexB1,gseEVexB2,gseEVexB3,gseEVexB4);
end
%% Assume normal electric fields are directly proportional to ion pressure gradient at boundary
c_eval('gseGradPi?_fromE = gseE?.resample(ne?)*ne?*units.e*1e-3*1e6*1e9*1e3; gseGradPi?_fromE.units = ''nPa/km'';',ic)
c_eval('Pi?perp = irf.ts_scalar(facPi?.time,(facPi?.yy.data+facPi?.zz.data)/2);',ic)
c_eval('Lp?= Pi?perp/ne?.resample(Pi?perp)/gseE?perp.abs/units.e*1e-9*1e-6*1e3;',ic)
c_eval('gseE?perp_filt = gseE?perp.filt(0,3,[],3);',ic)
c_eval('Lp?_filt= Pi?perp/ne?.resample(Pi?perp)/gseE?perp_filt.resample(Pi?perp).abs/units.e*1e-9*1e-6*1e3;',ic)

%%
disp('Done preparing data. Not MVA system.')