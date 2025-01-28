%% Prepare data
c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);
%mms_2015Oct16.construct_particle_data % Make omni fluxes
ic = 1:4;
avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
avNi = (ni1+ni2.resample(ni1.time)+ni3.resample(ni1.time)+ni4.resample(ni1.time))/4; avNi.name = '<ni>';
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; gseAvB.name = '<B> (GSE)';
gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; gseAvVe.name = '<ve> (GSE)';
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.name = '<J> (GSE)';
vDe = gseGradPe.cross(gseAvB.resample(gseGradPe))/units.e/avNe/gseAvB.abs2*1e-9*1e-3; vDe.units = 'km/s'; vDe.name = 'v_De';
vDi = -gseGradPe.cross(gseAvB.resample(gseGradPe))/units.e/avNi/gseAvB.abs2*1e-9*1e-3; vDi.units = 'km/s'; vDi.name = 'v_Di';

% Plasma beta and magnetic pressure
c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

% Non-ideal electric field
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

% Energy dissipation
c_eval('RedJ? = gseEVexB?.resample(gseJ?.time).dot(gseJ?)*1e-3; RedJ?.units = ''nW/m^3''; RedJ?.name = ''E+VexB d J'';',ic)
c_eval('Ge? = RedJ?*1e-9/(ne?*1e6)/fce?/units.e; Ge?.units = ''eV''; Ge?.name = ''Ge'';',ic)
c_eval('RedJe? = gseEVexB?.resample(gseJe?.time).dot(gseJe?)*1e-3; RedJe?.units = ''nW/m^3'';',ic)
c_eval('RedJi? = gseEVexB?.resample(gseJi?.time).dot(gseJi?)*1e-3; RedJi?.units = ''nW/m^3'';',ic)
c_eval('EdJ? = gseE?.resample(gseJ?.time).dot(gseJ?)*1e-3; EdJ?.units = ''nW/m^3'';',ic)
c_eval('EdJe? = gseE?.resample(gseJe?.time).dot(gseJe?)*1e-3; EdJe?.units = ''nW/m^3'';',ic)
c_eval('EdJi? = gseE?.resample(gseJi?.time).dot(gseJi?)*1e-3; EdJi?.units = ''nW/m^3'';',ic)
c_eval('EdVe? = gseE?.resample(gseVe?.time).dot(gseVe?); EdVe?.units = ''V/s'';',ic)
c_eval('gseEdJ?vec = gseE?.resample(gseJ?.time).*gseJ?*1e-3; EdJ?vec.units = ''nW/m^3'';',ic)
c_eval('gseRedJ?vec = gseEVexB?.resample(gseJ?.time).*(gseJ?)*1e-3; RedJ?vec.units = ''nW/m^3'';',ic)
c_eval('gseRedJe?vec = gseEVexB?.resample(gseJe?.time).*(gseJe?)*1e-3; RedJe?vec.units = ''nW/m^3'';',ic)
c_eval('gseRedJi?vec = gseEVexB?.resample(gseJi?.time).*(gseJi?)*1e-3; RedJi?vec.units = ''nW/m^3'';',ic)
c_eval('gseEdJ?vec = gseE?.resample(gseJ?.time).*gseJ?*1e-3; EdJ?vec.units = ''nW/m^3'';',ic)
c_eval('gseEdJe?vec = gseE?.resample(gseJe?.time).*(gseJe?)*1e-3; EdJe?vec.units = ''nW/m^3'';',ic)
c_eval('gseEdJi?vec = gseE?.resample(gseJi?.time).*(gseJi?)*1e-3; EdJi?vec.units = ''nW/m^3'';',ic)
c_eval('gseEdVe?vec = gseE?.resample(gseVe?.time).*(gseVe?); EdVe?vec.units = ''V/s'';',ic)
EJs = {'gseRedJ?vec','RedJe?vec','RedJi?vec','EdJ?vec','EdJe?vec','EdJi?vec'};

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% Calculate electric field from E dot B = 0
angle_lim = 20;
c_eval('[dslE?EdB0,Bangle?]=irf_edb(dslE?,dmpaB?,angle_lim,''Eperp+NaN'');',ic)
c_eval('gseE?EdB0 = mms_dsl2gse(dslE?EdB0,defatt?);',ic)

% Field aligned pressure and temperature tensors
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe? = irf.ts_tensor_xyz(facPe?.time,facPe?.data);',ic)
c_eval('facPe?pp = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp''); facPe?pp = irf.ts_tensor_xyz(facPe?pp.time,facPe?.data);',ic)
c_eval('facPe?qq = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq''); facPe?pp = irf.ts_tensor_xyz(facPe?qq.time,facPe?.data);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?); facTe? = irf.ts_tensor_xyz(facTe?.time,facTe?.data);',ic)

% Agyrotropy and anisotropy
c_eval('Q? = irf.ts_scalar(facPe?pp.time,(facPe?pp.data(:,1,2).^2+facPe?pp.data(:,1,3).^2+facPe?pp.data(:,2,3).^2)./(facPe?pp.data(:,2,2).^2+2*facPe?pp.data(:,2,2).*facPe?pp.data(:,1,1))); Q?.name = ''Q'';',ic);
c_eval('Dng? = irf.ts_scalar(facPe?pp.time,sqrt(8*(facPe?pp.data(:,1,2).^2+facPe?pp.data(:,1,3).^2+facPe?pp.data(:,2,3).^2))./(facPe?pp.data(:,1,1)+2*facPe?pp.data(:,2,2))); Dng?.name = ''Dng'';',ic);
% Compute agyrotropy Aphi from Pe?qq
c_eval('Ao? = irf.ts_scalar(facPe?.time,2*abs(facPe?qq.data(:,2,2)-facPe?qq.data(:,3,3))./(facPe?qq.data(:,2,2)+facPe?qq.data(:,3,3))); Ao?.name = ''Ao'';',ic);
% Compute temperature ratio An
c_eval('T?ratio = irf.ts_scalar(facTe?.time,facPe?pp.data(:,1,1)./(facPe?pp.data(:,2,2))); T?ratio.name = ''Tpar/Tperp'';',ic);
% Average work done on all the particles per gyroperiod per average thermal energy
c_eval('eEps? = abs(6*pi*EdVe?/(2*pi*fce?.resample(facTe?)*facTe?.trace)); eEps?.name = ''epsilon e'';',ic);

% Parallel and perpendicular velocity and field components
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)

c_eval('try [gseE?fastpar,gseE?fastperp] = irf_dec_parperp(gseB?.resample(gseE?fast),gseE?fast); end',ic); 
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?.resample(gseE?),gseE?);',ic)

c_eval('try gseVExB?fast = gseE?fast.cross(gseB?.resample(gseE?fast))/gseB?.resample(gseE?fast).abs2*1e3; gseVExB?fast.units = ''km/s''; end',ic);
c_eval('gseVExB? = gseE?.resample(gseB?).cross(gseB?)/gseB?.abs2*1e3; gseVExB?.units = ''km/s'';',ic)

% Prepare pitchangle 
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,13);',ic)

%% MVA: Rotate into lmn coordinates
ic=1:4;

avPe = (gsePe1.trace+gsePe2.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1))/3/4;
LgradP = avPe/gseGradPe.abs; LgradP.name = 'L_{P}';

% Proxy for magnetic field slippage (Scudder2105)
avFce = (fce1+fce2.resample(fce1)+fce3.resample(fce1)+fce4.resample(fce1))/4;
avRhoe = (re1+re2.resample(re1)+re3.resample(re1)+re4.resample(re1))/4;
Y = avFce*(avRhoe/LgradP).^2; Y.units = 's^-1';
Ynodim = (avRhoe/LgradP); Y.units = 's^-1';

%
facAvTe = (facTe1+facTe2.resample(facTe1)+facTe3.resample(facTe1)+facTe4.resample(facTe1))/4;
% Rotation of E+VexB (Scudder2105)
gseRotRe = mms_2015Oct16.rotRe(gseR1,gseR2,gseR3,gseR4,gseEVexB1,gseEVexB2,gseEVexB3,gseEVexB4);

% Set up coordinate system
% Default
[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
L = v(1,:); M = v(2,:); N = v(3,:);
coordLabels = {'L','M','N'};
lmn = [N;-M;L];

coordSystem = 1;
switch coordSystem % Choose another
  case 1 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 2 % N: minimum variance of J
    [out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-10-16T10:33:23.000Z/2015-10-16T10:33:32.000Z')));
    L = -v(2,:); M = -v(1,:); N = v(3,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
    [out,l,v] = irf_minvar(gseJcurl.tlim(irf.tint('2015-10-16T10:33:22.595Z/2015-10-16T10:33:31.284Z')));
    L = -v(1,:); M = v(2,:); N = -v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [L;M;N];        
  case 3 % N: maximum variance of E
    tint = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.400Z');
    [out,l,v] = irf_minvar(gseE3.tlim(tint));
    L = v(2,:); M = -v(3,:); N = -v(1,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
  case 4 % N: magnetosheath side normal derived from mms1 and mms4
    gseVec14 = gseR4-gseR1; gseVec14 = gseVec14.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec14.time,M);
    gseNorm14 = gseVec14.cross(gseM);
    gseNormalMSH = gseNorm14/gseNorm14.abs;

    N = -gseNormalMSH.data;
    M = M;
    L = cross(M,N);
  case 5 % N: magnetosphere side normal derived from mms3 and mms4
    gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec34.time,M);
    gseNorm34 = gseVec34.cross(gseM);
    gseNormalMSP = gseNorm34/gseNorm34.abs;

    N = -gseNormalMSP.data;
    M = M;
    L = cross(M,N);
end

% Rotate data
c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')
c_eval('mvaB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(L).data gseB?.dot(M).data gseB?.dot(N).data]); mvaB?.name = ''B LMN'';')
c_eval('mvaE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(L).data gseE?.dot(M).data gseE?.dot(N).data]); mvaE?.name = ''E LMN'';')
c_eval('mvaVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(L).data gseVe?.dot(M).data gseVe?.dot(N).data]); mvaVe?.name = ''Ve LMN'';')
c_eval('mvaVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(L).data gseVi?.dot(M).data gseVi?.dot(N).data]); mvaVi?.name = ''Vi LMN'';')
c_eval('mvaJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(L).data gseJ?.dot(M).data gseJ?.dot(N).data]); mvaJ?.units = gseJ?.units; mvaJ?.name = ''J LMN'';')
c_eval('mvaJe? = irf.ts_vec_xyz(gseJe?.time,[gseJe?.dot(L).data gseJe?.dot(M).data gseJe?.dot(N).data]); mvaJe?.units = gseJe?.units; mvaJe?.name = ''Je LMN'';')
c_eval('mvaJi? = irf.ts_vec_xyz(gseJi?.time,[gseJi?.dot(L).data gseJi?.dot(M).data gseJi?.dot(N).data]); mvaJi?.units = gseJi?.units; mvaJi?.name = ''Ji LMN'';')
mvaJcurl = irf.ts_vec_xyz(gseJcurl.time,[gseJcurl.dot(L).data gseJcurl.dot(M).data gseJcurl.dot(N).data]); mvaJcurl.units = gseJcurl.units;
c_eval('mvaPe? = mms.rotate_tensor(gsePe?,''rot'',L,M,N); mvaPe? = irf.ts_tensor_xyz(mvaPe?.time,mvaPe?.data); mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTe? = mms.rotate_tensor(gseTe?,''rot'',L,M,N); mvaTe? = irf.ts_tensor_xyz(mvaTe?.time,mvaTe?.data); mvaTe?.units = gseTe?.units;',ic)
c_eval('mvaVexB? = irf.ts_vec_xyz(gseVexB?.time,[gseVexB?.dot(L).data gseVexB?.dot(M).data gseVexB?.dot(N).data]); mvaVexB?.units = ''mV/m'';')
c_eval('mvaVixB? =  irf.ts_vec_xyz(gseVixB?.time,[gseVixB?.dot(L).data gseVixB?.dot(M).data gseVixB?.dot(N).data]); mvaVixB?.units = ''mV/m'';')
c_eval('mvaEVexB? =  irf.ts_vec_xyz(gseEVexB?.time,[gseEVexB?.dot(L).data gseEVexB?.dot(M).data gseEVexB?.dot(N).data]); mvaEVexB?.units = ''mV/m'';')
mvaVDe =  irf.ts_vec_xyz(vDe.time,[vDe.dot(L).data vDe.dot(M).data vDe.dot(N).data]); mvaVDe.units = '';
mvaAvJ =  irf.ts_vec_xyz(gseAvJ.time,[gseAvJ.dot(L).data gseAvJ.dot(M).data gseAvJ.dot(N).data]); mvaAvJ.units = 'nA/m^2';
c_eval('mvaJxB? = mvaJ?.cross(mvaB?.resample(mvaJ?.time));')
c_eval('mvaVExB? =  irf.ts_vec_xyz(gseVExB?.time,[gseVExB?.dot(L).data gseVExB?.dot(M).data gseVExB?.dot(N).data]);')
c_eval('mvaVe?par = gseVe?par;')
c_eval('mvaVe?perp = irf.ts_vec_xyz(gseVe?perp.time,[gseVe?perp.dot(L).data gseVe?perp.dot(M).data gseVe?perp.dot(N).data]);')
c_eval('mvaE?par = gseE?par;')
c_eval('mvaE?perp = irf.ts_vec_xyz(gseE?perp.time,[gseE?perp.dot(L).data gseE?perp.dot(M).data gseE?perp.dot(N).data]);')
%c_eval('mvaE?fastpar = gseE?fastpar;')
%c_eval('mvaE?fastperp = irf.ts_vec_xyz(gseE?fastperp.time,[gseE?fastperp.dot(L).data gseE?fastperp.dot(M).data gseE?fastperp.dot(N).data]);')
c_eval('mvaVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(L).data gseVi?.dot(M).data gseVi?.dot(N).data]);')

c_eval('mvaEdJ?vec = irf.ts_vec_xyz(gseEdJ?vec.time,[gseEdJ?vec.dot(L).data gseEdJ?vec.dot(M).data gseEdJ?vec.dot(N).data]);')
c_eval('mvaEdJe?vec = irf.ts_vec_xyz(gseEdJe?vec.time,[gseEdJe?vec.dot(L).data gseEdJe?vec.dot(M).data gseEdJe?vec.dot(N).data]);')
c_eval('mvaEdJi?vec = irf.ts_vec_xyz(gseEdJi?vec.time,[gseEdJi?vec.dot(L).data gseEdJi?vec.dot(M).data gseEdJi?vec.dot(N).data]);')
c_eval('mvaRedJ?vec = irf.ts_vec_xyz(gseRedJ?vec.time,[gseRedJ?vec.dot(L).data gseRedJ?vec.dot(M).data gseRedJ?vec.dot(N).data]);')
c_eval('mvaRedJe?vec = irf.ts_vec_xyz(gseRedJe?vec.time,[gseRedJe?vec.dot(L).data gseRedJe?vec.dot(M).data gseRedJe?vec.dot(N).data]);')
c_eval('mvaRedJi?vec = irf.ts_vec_xyz(gseRedJi?vec.time,[gseRedJi?vec.dot(L).data gseRedJi?vec.dot(M).data gseRedJi?vec.dot(N).data]);')

mvaRotRe = irf.ts_vec_xyz(gseRotRe.time,[gseRotRe.dot(L).data gseRotRe.dot(M).data gseRotRe.dot(N).data]);
mvaGradPe = irf.ts_vec_xyz(gseGradPe.time,[gseGradPe.dot(L).data gseGradPe.dot(M).data gseGradPe.dot(N).data]);

mvaAvE = (mvaE1+mvaE2.resample(mvaE1.time)+mvaE3.resample(mvaE1.time)+mvaE4.resample(mvaE1.time))/4; 
mvaAvVe = (mvaVe1+mvaVe2.resample(mvaVe1.time)+mvaVe3.resample(mvaVe1.time)+mvaVe4.resample(mvaVe1.time))/4; 
mvaAvVeperp = (mvaVe1perp+mvaVe2perp.resample(mvaVe1perp.time)+mvaVe3perp.resample(mvaVe1perp.time)+mvaVe4perp.resample(mvaVe1perp.time))/4; 
mvaAvVi = (mvaVi1+mvaVi2.resample(mvaVi1.time)+mvaVi3.resample(mvaVi1.time)+mvaVi4.resample(mvaVi1.time))/4; 
mvaAvB = (mvaB1+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4; 
mvaAvJ = (mvaJ1+mvaJ2.resample(mvaJ1.time)+mvaJ3.resample(mvaJ1.time)+mvaJ4.resample(mvaJ1.time))/4; mvaAvJ.units = mvaJ1.units;
mvaAvVExB = (mvaVExB1+mvaVExB2.resample(mvaVExB1.time)+mvaVExB3.resample(mvaVExB1.time)+mvaVExB4.resample(mvaVExB1.time))/4; 
mvaAvVexB = (mvaVexB1+mvaVexB2.resample(mvaVexB1.time)+mvaVexB3.resample(mvaVexB1.time)+mvaVexB4.resample(mvaVexB1.time))/4; 
mvaAvVixB = (mvaVixB1+mvaVixB2.resample(mvaVixB1.time)+mvaVixB3.resample(mvaVixB1.time)+mvaVixB4.resample(mvaVixB1.time))/4; 
mvaAvB = (mvaB1+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4;

mvaR0 = (mvaR1.resample(mvaR1.time)+mvaR2.resample(mvaR1.time)+mvaR3.resample(mvaR1.time)+mvaR4.resample(mvaR1.time))/4;
c_eval('mvaRR? = mvaR?-mvaR0; mvaRR? = mvaRR?.resample(irf_time(''2015-10-16T10:33:30.000Z'',''utc>epochTT'')).data;',ic)
c_eval('[mvaVe?par,mvaVe?perp,mvaVe?PA]=irf_dec_parperp(mvaB?.resample(mvaVe?),mvaVe?);',ic)
c_eval('[mvaVi?par,mvaVi?perp,mvaVi?PA]=irf_dec_parperp(mvaB?.resample(mvaVi?),mvaVi?);',ic)
 
%% Make Ohm's law terms
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

%% Figure 1 v1: Overview
tint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere

h = irf_plot(10);
ic = 1;
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.e64.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,18).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
end
if 1 % Ve
  c_eval('Veplot = gseVe?;',ic)
  c_eval('neplot = ne?;',ic)
  tmpdata = Veplot.data;
  indLowNe = find(neplot.data<1);
  tmpdata(indLowNe,:) = repmat([NaN NaN NaN],numel(indLowNe),1);
  Veplot.data =tmpdata;
  c_eval('Veplot = gseVe?.clone(gseVe?.time,tmpdata);',ic)
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{Veplot.x.tlim(tint),Veplot.y.tlim(tint),Veplot.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'n>1 cm^{-3}'},[0.08 0.95],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [-1200 700];
  
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hout,irf_colormap('space'))  
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.pitchangles(dmpaB?,24).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 0 % i DEF omni 64
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hout,irf_colormap('space'))  
end
if 0 % i DEF omni 32
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end

%load('caa/cmap.mat');
for ii = [4 5 7 8]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:10
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
end

irf_zoom(h,'x',tint)
irf_zoom(h([1:3 6 10]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure 1 v2: Overview: Ions larges scale, zoom in on electrons 
iTint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z');
eTint = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:35.00Z'); 
npanels = 13;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni error');  
  c_eval('irf_spectrogram(hca,iPDistError?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni counts');  
  c_eval('irf_spectrogram(hca,iPDistCounts?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  %hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 1 % iPDist pa 32
  iisub = iisub+1;
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,18).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap)
end

hca = irf_panel('delete for space');

if 0 % Ve n>1cc
  c_eval('Veplot = gseVe?;',ic)
  c_eval('neplot = ne?;',ic)
  tmpdata = Veplot.data;
  indLowNe = find(neplot.data<1);
  tmpdata(indLowNe,:) = repmat([NaN NaN NaN],numel(indLowNe),1);
  Veplot.data =tmpdata;
  c_eval('Veplot = gseVe?.clone(gseVe?.time,tmpdata);',ic)
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{Veplot.x.tlim(tint),Veplot.y.tlim(tint),Veplot.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'n>1 cm^{-3}'},[0.08 0.95],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [-1200 700];  
end
if 1 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.15],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux');  
  eint = [0 200];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [200 1000];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all e');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).pitchangles(dmpaB?,20).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

if 0 % i DEF omni 64
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hout,irf_colormap('space'))  
end
if 0 % i DEF omni 32
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.98 0.9],'fontsize',12,'color','k');
end

%load('caa/cmap.mat');
for ii = [4 5 8 9]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',eTint)
%irf_zoom(h([1:3 6 10]),'y')
hca = irf_panel('Te'); irf_zoom(hca,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
hca = irf_panel('delete for space'); delete(hca) 

%% Figure 1: Stripped version: LMN
iTint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z');
eTint = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z'); 
boundaryTint = EpochTT(['2015-10-16T10:33:26.20Z';...
                        '2015-10-16T10:33:27.00Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);
boundaryTint = EpochTT(['2015-10-16T10:33:27.30Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);                
npanels = 12;
cmap = 'jet';
h = irf_plot(npanels);
ic = 4;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?},''comp'');',ic)
  hca.YLabel.String = {'n_i','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 1 % eDist omni 64
  iisub = iisub+1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
  hhleg.FontSize = 9;
end

hca = irf_panel('delete for space');

if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  labelsOutside = 1;
  labelFontSize = 10;
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvE.z.resample(mvaOhmVexB.time),1*mvaOhmVexB.z,mvaOhmGradPe.z,mvaOhmVexB.z+mvaOhmGradPe.z.resample(mvaOhmVexB.time)+mvaAvE.z.resample(mvaOhmVexB.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('1')); hl = irf_legend(hca,{'E+v_{e}xB +\nabla\cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('1')); hl(1).VerticalAlignment = 'top';
    end
  else
    irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','v_{e}xB+\nabla \cdot P_e/ne+E'},[0.98 0.1],'fontsize',12);
  end
  
  %irf_legend(hca,{'4 sc average'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc average',ic),[0.06 0.9],'color','k','fontsize',11);    
end
if 1 % sqrt(Q)
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{sqrt(Q1),sqrt(Q2),sqrt(Q3),sqrt(Q4)},'comp');      
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,{'$$\sqrt{Q}$$',''},'interpreter','latex');
  %ylabel(hca,{'$$Q^{1/2}$$',''},'interpreter','latex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.06 0.9],'fontsize',11);
end

%load('caa/cmap.mat');
for ii = [4 5]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  h(ii).YLabel.FontSize = 12;
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',eTint)
%irf_zoom(h([1:3 6 10]),'y')


irf_pl_mark(h(1:iisub),eTint.epochUnix')
%irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix')

hca = irf_panel('Te'); irf_zoom(hca,'y');

hca = irf_panel('Vi'); hca.YLim = [-300 100];
hca = irf_panel('B'); hca.YLim = [-35 50];
hca = irf_panel('n'); hca.YLim = [0 29.9];
hca = irf_panel('J mom'); hca.YLim = [-1600 1200];
hca = irf_panel('Ve'); hca.YLim = [-800 700];
hca = irf_panel('Q'); hca.YLim = [0 0.09]; hca.YTick = [0 0.05];[0.04 0.08];
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-6 6];
hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'};
hca.XLabel.String = '';

ihmark = irf_pl_mark(h(1:iisub),eTint.epochUnix',[0.8 0.8 0.8]); 

hmark = irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark), hmark(ii).LineStyle = '-'; end

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

hca = irf_panel('delete for space'); 
hca.Visible = 'off';

% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(iisub+2);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.83 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';

%% Figure 1: Stripped version: LMN, average fields
iTint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z');
eTint = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z'); 
boundaryTint = EpochTT(['2015-10-16T10:33:26.20Z';...
                        '2015-10-16T10:33:27.00Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);
boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);                
npanels = 12;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{mvaAvB.x,mvaAvB.y,mvaAvB.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,{avNi,avNi},'comp');
  hca.YLabel.String = {'n_i','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('11'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_plot(hca,{mvaAvVi.x,mvaAvVi.y,mvaAvVi.z},'comp');
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
  hhleg=irf_legend(hca,irf_ssub('MMS ?',ic),[0.05 0.9],'color',[0 0 0],'fontsize',11);
  hcb.YLabel.String = {''};
end
if 1 % eDist omni 64
  iisub = iisub+1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
  hhleg.FontSize = 9;
  hhleg=irf_legend(hca,irf_ssub('MMS ?',ic),[0.05 0.9],'color',[0 0 0],'fontsize',11);
  hcb.YLabel.String = {'Differential Energy Flux',hcb.YLabel.String{2}};
  hcb.YLabel.Position(2)=hcb.YLabel.Position(2)+0.5;
  hcb.YLabel.FontSize = 9;
end

hca = irf_panel('delete for space');

if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{mvaAvJ.x.tlim(tint),mvaAvJ.y.tlim(tint),mvaAvJ.z.tlim(tint)},'comp');
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{mvaAvVe.x.tlim(tint),mvaAvVe.y.tlim(tint),mvaAvVe.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  irf_plot(hca,{facAvTe.xx.tlim(tint),(facAvTe.yy+facAvTe.zz)/2},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  labelsOutside = 1;
  labelFontSize = 10;
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvE.z.resample(mvaOhmVexB.time),-1*mvaOhmVexB.z,mvaOhmGradPe.z,mvaOhmVexB.z+mvaOhmGradPe.z.resample(mvaOhmVexB.time)+mvaAvE.z.resample(mvaOhmVexB.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','-v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('1')); hl = irf_legend(hca,{'E+v_{e}xB +\nabla\cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('1')); hl(1).VerticalAlignment = 'top';
    end
  else
    irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','v_{e}xB+\nabla \cdot P_e/ne+E'},[0.98 0.1],'fontsize',12);
  end
  
  %irf_legend(hca,{'4 sc average'},[0.98 0.1],'fontsize',12);  
end
if 1 % sqrt(Q)
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{sqrt(Q1),sqrt(Q2),sqrt(Q3),sqrt(Q4)},'comp');      
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,{'$$\sqrt{Q}$$',''},'interpreter','latex');
  %ylabel(hca,{'$$Q^{1/2}$$',''},'interpreter','latex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.06 0.9],'fontsize',11);
end

%load('caa/cmap.mat');
for ii = [4 5]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',11)
  h(ii).YLabel.FontSize = 11;
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',eTint)
%irf_zoom(h([1:3 6 10]),'y')


irf_pl_mark(h(1:iisub),eTint.epochUnix')
%irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix')

hca = irf_panel('Te'); irf_zoom(hca,'y');

hca = irf_panel('Vi'); hca.YLim = [-300 100];
hca = irf_panel('B'); hca.YLim = [-35 50];
hca = irf_panel('n'); hca.YLim = [0 29.9];
hca = irf_panel('J mom'); hca.YLim = [-1600 1200];
hca = irf_panel('Ve'); hca.YLim = [-800 700];
hca = irf_panel('Q'); hca.YLim = [0 0.09]; hca.YTick = [0 0.05];[0.04 0.08];
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-6 6];
hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'};
hca.XLabel.String = '';

ihmark = irf_pl_mark(h(1:iisub),eTint.epochUnix',[0.8 0.8 0.8]); 

% Mark magnetospheric separatrix
hmark1 = irf_pl_mark(h(iisub+2:npanels),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; end

% Mark magnetospheath separatrix
hmark2 = irf_pl_mark(h(iisub+2:npanels),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; end


irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

hca = irf_panel('delete for space'); 
hca.Visible = 'off';

% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(iisub+2);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center'; hleg_mspsep.FontSize = 11;
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center'; hleg_outflow.FontSize = 11;
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center'; hleg_mshsep.FontSize = 11;

%% Figure 2: Plot figure with fields etc of 4 sc, for paper
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
npanels = 10;

boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);                
pshift = 3;                      
h = irf_plot(npanels + pshift);
irf_panel('delete1');
irf_panel('delete2');
irf_panel('delete3');
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{mvaEVexB1.z,mvaEVexB2.z,mvaEVexB3.z,mvaEVexB4.z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  hca.YLim = [-10 10];
end
if 1 % (E + vexB)*J, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  %c_eval('filtEdJ? = RedJ?;'); 
  hca = irf_panel('E'' dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{RedJ1,RedJ2,RedJ3,RedJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  %hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  %hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('E dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJ1,filtEdJ2,filtEdJ3,filtEdJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = RedJe?;')
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP    
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h(4:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

drawnow

% reset ylims 
hca = irf_panel('BN'); hca.YLim = [-7 7];
hca = irf_panel('BM'); hca.YTick = [-10 0 10];
hca = irf_panel('JM'); hca.YTick = [-1000 0 1000];
hca = irf_panel('E + vexB'); hca.YLim = [-4.5 6.5];
hca = irf_panel('EN'); hca.YLim = [-7.5 6.5];


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
shift = 3; % three deleted panels
legshift = 2; % the two sc configuration plots
for ii = 1:npanels
  irf_legend(h(ii+shift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+shift).FontSize = 12;  
  %if labelsOutside
  %  h(ii+shift).Position(3) = h(ii+shift).Position(3)*0.88;
  %end  
  h(ii+shift).YLabel.FontSize = 11;
end

for ii = 1:3; h(ii).Visible = 'off'; end

% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(pshift+1);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';

% Add borders/separtrices for the different regions
hmark1 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
hmark2 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; hmark1(ii).Color = [0.4 0.4 0.4]; end
for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; hmark2(ii).Color = [0.4 0.4 0.4]; end




%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

% Plot sc positions
if exist('h2'); delete(h2); end
nrows = 5;
ncols = 3;
h2(1) = subplot(nrows,ncols,1); %h2(1).Position(2) = h2(1).Position(2)+0.02;
h2(2) = subplot(nrows,ncols,2); %h2(2).Position(2) = h2(2).Position(2)+0.02;

mms_marker={{'ks','markersize',10},{'rd','markersize',10},...
	{'go','markersize',10,'color',[0 0.6 0]},{'bv','markersize',10}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
	{'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
	{'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
sc_list = 1:4;

x = {mvaRR1,mvaRR2,mvaRR3,mvaRR4};

hca = h2(1);
hold(hca,'on');
for ic=1:4
  % put Cluster markers  
  plot(hca,x{ic}(3),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.5*[-10 10];
hca.YLim = 1.5*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'L (km)';

%plot(mvaRR1(3),mvaRR1(1))
hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.7*[-10 10];
hca.YLim = 1.7*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'M (km)';
hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
hleg.Box = 'off';
hleg.Position(1) = hleg.Position(1)+0.2;
hleg.Position(2) = hleg.Position(2)+0.05;

% plot plane
v = [-0.90 0.26 0.36];
v = [-0.90 -0.28 -0.33];
lmnV = [-0.57 -0.03 -0.82];
lmnV = [-0.55 -0.08 -0.83];
lmnV = [-0.45 -0.05 -0.89];
lmnV = [-0.49 -0.05 -0.87];
v = [-0.9139   -0.4006    0.0653];
v = [-0.90 -0.28 -0.33];

mshV = 55*[-0.88 -0.26 -0.40]; % GSE

mspV= 31.8*[-0.88 -0.42 0.24]; % GSE
mshV = 55*[-0.90 -0.28 -0.33]; %GSE
mspVlmn = mspV*[L' M' N'];
mshVlmn = mshV*[L' M' N'];
mspN = irf_norm(mspVlmn);
mshN = irf_norm(mshVlmn);

x = 50*[-1 1];
y = 50*[-1 1];
z = 50*[-1 1];

%zFun = @(x,y,n) -(n(1)*x+n(2)*y)/n(3);
% (planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3) = 0;
%  planeNormal(1)*x+planeNormal(2)*y = 0;
% x = - planeNormal(2)*y/planeNormal(1);
funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);



if exist('hmshN1'); delete(hmshN1); end
if exist('hmshN2'); delete(hmshN2); end
if exist('hmspN1'); delete(hmspN1); end
if exist('hmspN2'); delete(hmspN2); end

hca = h2(1); subplot(nrows,ncols,1);
hold(hca,'on');
%hmshN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN1 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
hmshN1 = plot(hca,funZ(x,y,mshN),x+18,'k-.');
hmspN1 = plot(hca,funZ(x,y,mspN)-10,x,'k-');

ht = text(9.5,9,'MSH'); ht.Rotation = 55; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
ht = text(-11.2,8,'MSP'); ht.Rotation = -80; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;

hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
%hmspN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN2 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
%hmspN1 = plot(hca,funZ(x,y,mshN),funY(x,z,mshN),'r');
%hmspN2 = plot(hca,funZ(x,y,mshN),funY(x,z,mspN),'b');
hmspN1 = plot(hca,z+10,funY(x,z,mshN),'k-.');
hmspN2 = plot(hca,z-10,funY(x,z,mspN),'k-');
hold(hca,'off');
%h2(1).Position(2) = h2(1).Position(2)+0.02;
%h2(2).Position(2) = h2(2).Position(2)+0.02;
irf_legend(h2(1),'a)',[-0.4 1],'color',[0 0 0])
irf_legend(h2(2),'b)',[-0.4 1],'color',[0 0 0])
%
for ii = 1:2;
  h2(ii).YLabel.FontSize = 12;
  h2(ii).XLabel.FontSize = 12;
  h2(ii).FontSize = 12;
  h2(ii).Position(2) = h2(ii).Position(2)+0.05;
  h2(ii).Position(1) = h2(ii).Position(1)+0.05;
end

%% Figure 2: Plot figure with fields etc of 4 sc
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
npanels = 12;
h = irf_plot(npanels + 3);
irf_panel('delete1');
irf_panel('delete2');
irf_panel('delete3');
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{mvaEVexB1.z,mvaEVexB2.z,mvaEVexB3.z,mvaEVexB4.z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  hca.YLim = [-10 10];
end

if 1 % (E + vexB)*J, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  %c_eval('filtEdJ? = RedJ?;'); 
  hca = irf_panel('E'' dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{RedJ1,RedJ2,RedJ3,RedJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  %hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  %hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('E dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJ1,filtEdJ2,filtEdJ3,filtEdJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = RedJe?;')
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP    
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 1 % sqrt(Q)
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{sqrt(Q1),sqrt(Q2),sqrt(Q3),sqrt(Q4)},'comp');    
  hca.YLabel.String = {'\sqrt{Q}',''};
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,{'$$\sqrt{Q}$$',''},'interpreter','latex');
  %ylabel(hca,{'$$Q^{1/2}$$',''},'interpreter','latex');
end
if 0 % Dng
  hca = irf_panel('Dng');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Dng1.tlim(tint),Dng2.tlim(tint),Dng3.tlim(tint),Dng4.tlim(tint)},'comp');    
  hca.YLabel.String = {'D_{ng,e}',''};
  %ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  labelsOutside = 1;
  labelFontSize = 10;
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_plot(hca,{mvaAvE.z,-1*mvaOhmVexB.z,-1*mvaOhmGradPe.z,-1*mvaOhmVexB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','-v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('b')); hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('b')); hl(1).VerticalAlignment = 'top';
    end
  else
    irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  end
  
  %irf_legend(hca,{'4 sc average'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc average',ic),[0.98 0.95],'color','k','fontsize',10);
  hca.YLim = [-10 10];
end

if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h(4:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

drawnow

% reset ylims 
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-5.5 5.5];
%hca = irf_panel('Dng'); hca.YLim = [0 0.15];
hca = irf_panel('Q'); hca.YLim = [0 0.09]; hca.YTick = [ 0.04 0.08];
hca = irf_panel('E + vexB'); hca.YLim = [-4.5 6.5];
hca = irf_panel('EN'); hca.YLim = [-7.5 6.5];
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-5 5];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
shift = 3; % three deleted panels
legshift = 2; % the two sc configuration plots
for ii = 1:npanels
  irf_legend(h(ii+shift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+shift).FontSize = 12;  
  if labelsOutside
    h(ii+shift).Position(3) = h(ii+shift).Position(3)*0.88;
  end  
  h(ii+shift).YLabel.FontSize = 11;
end

for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Plot sc positions
if exist('h2'); delete(h2); end
nrows = 5;
ncols = 3;
h2(1) = subplot(nrows,ncols,1); %h2(1).Position(2) = h2(1).Position(2)+0.02;
h2(2) = subplot(nrows,ncols,2); %h2(2).Position(2) = h2(2).Position(2)+0.02;

mms_marker={{'ks','markersize',10},{'rd','markersize',10},...
	{'go','markersize',10,'color',[0 0.6 0]},{'bv','markersize',10}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
	{'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
	{'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
sc_list = 1:4;

x = {mvaRR1,mvaRR2,mvaRR3,mvaRR4};

hca = h2(1);
hold(hca,'on');
for ic=1:4
  % put Cluster markers  
  plot(hca,x{ic}(3),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.5*[-10 10];
hca.YLim = 1.5*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'L (km)';

%plot(mvaRR1(3),mvaRR1(1))
hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.7*[-10 10];
hca.YLim = 1.7*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N';
hca.YLabel.String = 'M';
hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
hleg.Box = 'off';
hleg.Position(1) = hleg.Position(1)+0.2;

% plot plane
v = [-0.90 0.26 0.36];
v = [-0.90 -0.28 -0.33];
lmnV = [-0.57 -0.03 -0.82];
lmnV = [-0.55 -0.08 -0.83];
lmnV = [-0.45 -0.05 -0.89];
lmnV = [-0.49 -0.05 -0.87];
v = [-0.9139   -0.4006    0.0653];
v = [-0.90 -0.28 -0.33];

mshV = 55*[-0.88 -0.26 -0.40]; % GSE

mspV= 31.8*[-0.88 -0.42 0.24]; % GSE
mshV = 55*[-0.90 -0.28 -0.33]; %GSE
mspVlmn = mspV*[L' M' N'];
mshVlmn = mshV*[L' M' N'];
mspN = irf_norm(mspVlmn);
mshN = irf_norm(mshVlmn);

x = 50*[-1 1];
y = 50*[-1 1];
z = 50*[-1 1];

%zFun = @(x,y,n) -(n(1)*x+n(2)*y)/n(3);
% (planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3) = 0;
%  planeNormal(1)*x+planeNormal(2)*y = 0;
% x = - planeNormal(2)*y/planeNormal(1);
funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);



if exist('hmshN1'); delete(hmshN1); end
if exist('hmshN2'); delete(hmshN2); end
if exist('hmspN1'); delete(hmspN1); end
if exist('hmspN2'); delete(hmspN2); end

hca = h2(1); subplot(nrows,ncols,1);
hold(hca,'on');
%hmshN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN1 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
hmshN1 = plot(hca,funZ(x,y,mshN),x+18,'k-.');
hmspN1 = plot(hca,funZ(x,y,mspN)-10,x,'k-');

ht = text(9.5,9,'MSH'); ht.Rotation = 55; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
ht = text(-11.2,8,'MSP'); ht.Rotation = -80; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;

hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
%hmspN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN2 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
%hmspN1 = plot(hca,funZ(x,y,mshN),funY(x,z,mshN),'r');
%hmspN2 = plot(hca,funZ(x,y,mshN),funY(x,z,mspN),'b');
hmspN1 = plot(hca,z+10,funY(x,z,mshN),'k-.');
hmspN2 = plot(hca,z-10,funY(x,z,mspN),'k-');
hold(hca,'off');
%h2(1).Position(2) = h2(1).Position(2)+0.02;
%h2(2).Position(2) = h2(2).Position(2)+0.02;
irf_legend(h2(1),'a)',[-0.4 1],'color',[0 0 0])
irf_legend(h2(2),'b)',[-0.4 1],'color',[0 0 0])
%
for ii = 1:2;
  h2(ii).YLabel.FontSize = 12;
  h2(ii).XLabel.FontSize = 12;
  h2(ii).FontSize = 12;
  h2(ii).Position(2) = h2(ii).Position(2)+0.03;
end

%% Figure 2: Different comparison plot of more field components
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(13);
%irf_panel('delete1');
%irf_panel('delete2');
%irf_panel('delete3');
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).x,mvaE2.tlim(tint).x,mvaE3.tlim(tint).x,mvaE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).y,mvaE2.tlim(tint).y,mvaE3.tlim(tint).y,mvaE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, L 4sc
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.x+1*mvaVexB1.resample(mvaE1.time).x,...
                mvaE2.x+1*mvaVexB2.resample(mvaE2.time).x,...
                mvaE3.x+1*mvaVexB3.resample(mvaE3.time).x,...
                mvaE4.x+1*mvaVexB4.resample(mvaE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};      
  hca.YLabel.String = {'R_e_L','(mV/m)'};      
end
if 1 % E + vexB, M 4sc
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.y+1*mvaVexB1.resample(mvaE1.time).y,...
                mvaE2.y+1*mvaVexB2.resample(mvaE2.time).y,...
                mvaE3.y+1*mvaVexB3.resample(mvaE3.time).y,...
                mvaE4.y+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  hca.YLabel.String = {'R_e_M','(mV/m)'};            
end
if 1 % E + vexB, N 4sc
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};   
  hca.YLabel.String = {'R_e_N','(mV/m)'};         
end
if 1 % JeL
  hca = irf_panel('JeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).x,mvaJe2.tlim(tint).x,mvaJe3.tlim(tint).x,mvaJe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_e_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeM
  hca = irf_panel('JeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).y,mvaJe2.tlim(tint).y,mvaJe3.tlim(tint).y,mvaJe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_e_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeN
  hca = irf_panel('JeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).z,mvaJe2.tlim(tint).z,mvaJe3.tlim(tint).z,mvaJe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_e_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);  
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = EdJe?;',ic)
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'R_e\cdot J_e','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP  
  
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(mvaJe)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:9
  irf_legend(h(ii+3),legends{ii},[0.01 0.9],'color',[0 0 0])
end

%delete(h(1:3))

%% Figure 3: Ohm's law etc
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(10);
labelsOutside = 1; labelFontSize = 12;
%irf_panel('delete1');
%irf_panel('delete2');
%irf_panel('delete3');
if 1 % average B
  hca = irf_panel('av B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_plot(hca,{mvaAvB.x.tlim(tint),mvaAvB.y.tlim(tint),mvaAvB.z.tlim(tint),mvaAvB.abs.tlim(tint)},'comp');
  irf_plot(hca,{mvaAvB.x.tlim(tint),mvaAvB.y.tlim(tint),mvaAvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % average E
  hca = irf_panel('av E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_plot(hca,{mvaAvE.x.tlim(tint),mvaAvE.y.tlim(tint),mvaAvE.z.tlim(tint),mvaAvE.abs.tlim(tint)},'comp');
  irf_plot(hca,{mvaAvE.x.tlim(tint),mvaAvE.y.tlim(tint),mvaAvE.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'E_L','E_M','E_N','|E|'},[0.98 0.9],'fontsize',12);
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % gradPe Ohm
  hca = irf_panel('-grad Pe/ne');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{-mvaOhmGradPe.x.tlim(tint),-mvaOhmGradPe.y.tlim(tint),-mvaOhmGradPe.z.tlim(tint)},'comp');
  hca.YLabel.String = {'-\nabla \cdot P_e/ne','(mV/m)'};
end
if 1 % Ve x B Ohm
  hca = irf_panel('-VexB 4sc av');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{-mvaOhmVexB.x.tlim(tint),-mvaOhmVexB.y.tlim(tint),-mvaOhmVexB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e \times B','(mV/m)'};
end
if 1 % Vi x B Ohm
  hca = irf_panel('-VixB 4sc av');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{-mvaOhmVixB.x.tlim(tint),-mvaOhmVixB.y.tlim(tint),-mvaOhmVixB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_i \times B','(mV/m)'};
end
if 1 % J x B Ohm
  hca = irf_panel('JxB 4sc av');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaOhmJxB.x.tlim(tint),mvaOhmJxB.y.tlim(tint),mvaOhmJxB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'J \times B/ne','(mV/m)'};
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_plot(hca,{mvaAvE.z,-1*mvaOhmVexB.z,-1*mvaOhmGradPe.z,-1*mvaOhmVexB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','-v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('b')); hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('b')); hl(1).VerticalAlignment = 'top';
    end
  else
    
    irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  end
  %irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);
  hca.YLim = [-10 10];
end
if 1 % Ohms law ions: N
  hca = irf_panel('Ohms law ions: N');
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_plot(hca,{mvaAvE.z,mvaOhmJxB.z,-1*mvaOhmVixB.z,-1*mvaOhmGradPe.z,mvaOhmJxB.z.resample(mvaAvE.time)-1*mvaOhmVixB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','jxB/ne'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('az')); hl = irf_legend(hca,{'-v_{i}xB','-\nabla P_e/ne'},[1.01 0.65]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      hl = irf_legend(hca,{'jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[1.01 0.3],'color',mms_colors('b')); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x'));
      hl = irf_legend(hca,{'jxB/ne'},[1.01 0.79],'color',mms_colors('y'));
      hl = irf_legend(hca,{'-v_{i}xB'},[1.01 0.59],'color',mms_colors('a'));
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('z'));
      hl = irf_legend(hca,{'jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[1.01 0.19],'color',mms_colors('b'));    
    end
  else
    set(hca,'ColorOrder',mms_colors('xyazb'))  
    irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','JxB/ne-v_{i}xB-\nabla P_e/ne'},[0.98 0.1],'fontsize',12);
  end   
  
  %irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);
  hca.YLim = [-10 10];
end
if 1 % Ohms law 2: N
  hca = irf_panel('Ohms law 2: N');
  set(hca,'ColorOrder',mms_colors('xz'))
  irf_plot(hca,{mvaAvE.z+1*mvaOhmVexB.resample(mvaAvE.time).z,-1*mvaOhmGradPe.z},'comp'); 
  hca.YLabel.String = {'E_N','(mV/m)'};
  if labelsOutside    
    hl = irf_legend(hca,{'E+v_{e}xB'},[1.01 0.85],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';   hl(1).FontSize = labelFontSize;
    hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.5],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
  else
    set(hca,'ColorOrder',mms_colors('xz'))  
    irf_legend(hca,{'E+v_{e}xB','-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  end  
  %irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);  
  hca.YLim = [-10 10];
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = RedJe?;',ic)
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'R_e\cdot J_e','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % E*J, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = EdJ?;',ic)
  hca = irf_panel('E dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 1 % Y: rhoe/LgradP    
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(mvaJe)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:10
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  if labelsOutside
    h(ii).Position(3) = h(ii).Position(3)*0.8;
  end  
end
for ii = 1:10
  h(ii).FontSize = 12;
end
for ii = 2:9;
  h(ii).YLim = 0.95*[-5 5];
end

%% Figure: Energy dissipation plot
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(14);
%irf_panel('delete1');
%irf_panel('delete2');
%irf_panel('delete3');
if 1 % |B|
  hca = irf_panel('|B|');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 0 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tintZoom).x,mvaE2.tlim(tintZoom).x,mvaE3.tlim(tintZoom).x,mvaE4.tlim(tintZoom).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tintZoom).y,mvaE2.tlim(tintZoom).y,mvaE3.tlim(tintZoom).y,mvaE4.tlim(tintZoom).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tintZoom).z,mvaE2.tlim(tintZoom).z,mvaE3.tlim(tintZoom).z,mvaE4.tlim(tintZoom).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, L 4sc
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
%   irf_plot(hca,{mvaE1.x+1*mvaVexB1.resample(mvaE1.time).x,...
%                 mvaE2.x+1*mvaVexB2.resample(mvaE2.time).x,...
%                 mvaE3.x+1*mvaVexB3.resample(mvaE3.time).x,...
%                 mvaE4.x+1*mvaVexB4.resample(mvaE4.time).x},'comp'); 
  irf_plot(hca,{mvaE1.x+1*mvaVexB1.resample(mvaE1.time).x,...
                mvaE2.x+1*mvaVexB2.resample(mvaE2.time).x,...
                mvaE3.x+1*mvaVexB3.resample(mvaE3.time).x,...
                mvaE4.x+1*mvaVexB4.resample(mvaE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};      
  hca.YLabel.String = {'E''_L','(mV/m)'};      
end
if 1 % E + vexB, M 4sc
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.y+1*mvaVexB1.resample(mvaE1.time).y,...
                mvaE2.y+1*mvaVexB2.resample(mvaE2.time).y,...
                mvaE3.y+1*mvaVexB3.resample(mvaE3.time).y,...
                mvaE4.y+1*mvaVexB4.resample(mvaE4.time).y},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  hca.YLabel.String = {'E''_M','(mV/m)'};            
end
if 1 % E + vexB, N 4sc
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};   
  hca.YLabel.String = {'E''_N','(mV/m)'};         
end
if 0 % JeL
  hca = irf_panel('JeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).x,mvaJe2.tlim(tint).x,mvaJe3.tlim(tint).x,mvaJe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_e_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % JeM
  hca = irf_panel('JeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).y,mvaJe2.tlim(tint).y,mvaJe3.tlim(tint).y,mvaJe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_e_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % JeN
  hca = irf_panel('JeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).z,mvaJe2.tlim(tint).z,mvaJe3.tlim(tint).z,mvaJe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_e_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JL
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJ1.tlim(tint).x,mvaJ2.tlim(tint).x,mvaJ3.tlim(tint).x,mvaJ4.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJ1.tlim(tint).y,mvaJ2.tlim(tint).y,mvaJ3.tlim(tint).y,mvaJ4.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JN
  hca = irf_panel('JN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJ1.tlim(tint).z,mvaJ2.tlim(tint).z,mvaJ3.tlim(tint).z,mvaJ4.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);  
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % E*J L
  hca = irf_panel('E dot J L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaEdJ1vec.x,mvaEdJ2vec.x,mvaEdJ3vec.x,mvaEdJ4vec.x},'comp'); 
  hca.YLabel.String = {'E_L J_L','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J M
  hca = irf_panel('E dot J M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaEdJ1vec.y,mvaEdJ2vec.y,mvaEdJ3vec.y,mvaEdJ4vec.y},'comp'); 
  hca.YLabel.String = {'E_M J_M','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J N
  hca = irf_panel('E dot J N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaEdJ1vec.z,mvaEdJ2vec.z,mvaEdJ3vec.z,mvaEdJ4vec.z},'comp'); 
  hca.YLabel.String = {'E_N J_N','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J 
  hca = irf_panel('E dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{EdJ1,EdJ2,EdJ3,EdJ4},'comp'); 
  hca.YLabel.String = {'|E\cdot J|','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J L
  hca = irf_panel('E'' dot J L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaRedJ1vec.x,mvaRedJ2vec.x,mvaRedJ3vec.x,mvaRedJ4vec.x},'comp'); 
  hca.YLabel.String = {'E_L J_L','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J M
  hca = irf_panel('E'' dot J M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaRedJ1vec.y,mvaRedJ2vec.y,mvaRedJ3vec.y,mvaRedJ4vec.y},'comp'); 
  hca.YLabel.String = {'E_M J_M','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J N
  hca = irf_panel('E'' dot J N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaRedJ1vec.z,mvaRedJ2vec.z,mvaRedJ3vec.z,mvaRedJ4vec.z},'comp'); 
  hca.YLabel.String = {'E_N J_N','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 1 % E*J 
  hca = irf_panel('E'' dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{RedJ1,RedJ2,RedJ3,RedJ4},'comp'); 
  hca.YLabel.String = {'|E\cdot J|','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234')) 
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP  
  
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(mvaJe)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:9
  irf_legend(h(ii+3),legends{ii},[0.01 0.9],'color',[0 0 0])
end

%delete(h(1:3))

%% Figure 4: For paper, 1 sc, maybe 2 sc, horizontal reversed time positioning
scList = [1 2];
tintZoom = irf.tint('2015-10-16T10:33:28.00Z/2015-10-16T10:33:31.00Z');

nsc = numel(scList);

% Initialize timeseries plot
h1 = irf_plot(nsc);
if nsc == 1
  h1.Position(2) = 0.85; h1.Position(4) = 0.08;
elseif nsc == 2
  h1(1).Position(2) = 0.88; h1(1).Position(4) = 0.06;
  h1(2).Position(2) = 0.82; h1(2).Position(4) = 0.06;
end
    

for ic = scList
  hca = irf_panel(irf_ssub('Ve?',ic));
  c_eval('irf_plot(hca,{mvaVe?par,mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)
  irf_legend(hca,{'v_{||}','v_{\perp,L}','v_{\perp,M}','v_{\perp,N}'},[0.02 0.02])
  hca.YLabel.String = {'v_e','(km/s)'};
  irf_legend(hca,{irf_ssub('MMS ?',ic)},[0.1 0.95])
end
irf_zoom(h1,'x',tintZoom)
irf_zoom(h1,'y')
irf_plot_axis_align

% Initialize particle distribution plot
nRows = 4;
nCols = 5;
isub = 1;
clear h2;
for ii = nCols+1:(nRows+1)*nCols; h2(isub) = subplot(nRows+1,nCols,ii); isub = isub + 1; end

% for ii = 2:nRows+1
%   for jj = 1:nCols
%     isub = isub + 1;         
%     h2(isub) = subplot(nRows+1,nCols,jj+(ii-1)*totCols);
%     %h2left = [h2left; h2(isub).Position(1)];
%   end
% end

if 0 % shift positions
xLeft = h2(1).Position(1); xRight = h2(4).Position(1)+h2(4).Position(3);
xGap = 0.00;
xWidth = (xRight-xLeft)/ncols - (ncols/2-1)*xGap;

yTop = h2(1).Position(2)+h2(1).Position(3); yBot = h2(end).Position(2);
yGap = 0.05;
yWidth = (yTop-yBot)/nrows - (nrows/2-1)*yGap;
xWidth = yWidth;
ii = 1;
for irow = 1:nrows
  for icol = 1:ncols          
    xgap = ceil(icol/2-1)*xGap;
    x0 = xLeft + xgap + (icol-1)*xWidth;
    ygap = ceil(irow/2-1)*yGap;
    y0 = yTop - yWidth - ygap - (irow-1)*yWidth;
        
    pos = [x0 y0 xWidth yWidth];
    h2(ii).Position = pos;
    ii = ii + 1;    
  end
end
end

% Decide times
indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
indTimes = {[1214 1212 1210 1208 1206]+4,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
isub = 1;
for iicc = 1:numel(scList)
  ic = scList(iicc);
  c_eval('dist = ePDist?;',ic)
  times = dist.time;
  indTime = indTimes{ic};
  % Projection coordinate system
  csys = 1;
  switch csys
    case 1 % z: B, x: N, y: zxx    
      c_eval('z = gseB?.resample(ePDist?);',ic); 
      z = -z/z.abs;
      tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
      tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

      x = cross(z,cross(tsN,z));    
      y1 = cross(z,x);
      y2 = cross(z,cross(tsL,z));    
      y = y2;
      vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
    case 2
      x = hatB0;
      y = hatExB0;
      z = cross(x,y);
      vlabels = {'B','E\times B','B\times(E\times B)'};
    case 3
      x = [1 0 0];
      y = [0 1 0];
      z = [0 0 1];
      vlabels = {'X','Y','Z'};
    case 4
      x = L;
      y = M;
      z = N;
      vlabels = {'L','M','N'};
  end
  X = x;
  Y = y;
  Z = z;

  c_eval('dist = ePDist?;',ic)
  vlim = 15*1e3;
  elevlim = 10;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  haveYLabel = 0;
  for ii = 1:numel(indTime) % plot mms1, plane: NxB, N 
    x = X(indTime(ii)).data;
    y = Y(indTime(ii)).data;
    z = Z(indTime(ii)).data;

    time = times(indTime(ii));
    timeUTC = time.utc;  
    hmark = irf_pl_mark(h1(iicc),[time.epochUnix],mms_colors(irf_ssub('?',ic)));

    % Get mean vectors
    c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
    hatVe0 = double(irf_norm(Ve0));    
    c_eval('B0 = gseB?.resample(time).data;',ic); 
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = gseE?.resample(time).data;',ic); 
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);
   % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};


    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);    
    %hca.Title.String = timeUTC(15:23);
    hca.Title.String = '';
    ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    ht.VerticalAlignment = 'top';
    ht.HorizontalAlignment = 'left';
    ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    hca.XAxisLocation = 'top';


    if ~haveYLabel          
      hca.YLabel.String = 'N_{\perp}\times B';    
      haveYLabel = 1;
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end
    %hca.XTickLabel = ''; 
    %hca.XLabel.String = ''; 
    hca.XLabel.String = 'N_{\perp}';
    if 0
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[-y;-z;-x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
      titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
      hca.Title.String = titleStr;
      hca.Title.String = '';
      colormap(hca,strCMap)
      hca.XDir = 'reverse';

      hca.XLabel.String = 'N_{\perp}';
      hca.YLabel.String = 'B';
    end
  end 
  haveYLabel = 0;
  for ii = 1:numel(indTime) % plot mms1, plane: B, ... 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;
    
  time = times(indTime(ii));
  timeUTC = time.utc;  
  hmark = irf_pl_mark(h1(iicc),[time.epochUnix],mms_colors(irf_ssub('?',ic)));
  
  % Get mean vectors
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

  
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[x;-z;-y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
  mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  %hca.XDir = 'reverse';

  hca.XLabel.String = 'N_{\perp}\times B';
  hca.YLabel.String = 'B';
  
  
  if ~haveYLabel          
    hca.YLabel.String = 'B';    
    haveYLabel = 1;
  else
    hca.YLabel.String = '';
    hca.YTickLabel = ''; 
  end
  %hca.XTickLabel = ''; 
  %hca.XLabel.String = ''; 
  %hca.XLabel.String = 'N_{\perp}';

end 
end
for ii = 1:nRows*nCols
  originalPosition{ii} = h2(ii).Position;
end

hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
%hCB.Position = [h2(nCols).Position(1)+xWidth+0.02 h2(4).Position(2) 0.02 h2(4).Position(4)];

cmap = irf_colormap('space');
cmap = 'jet';
for ii = 1:nRows*nCols
  h2(ii).Position(3) = originalPosition{ii}(3)*2.7;
  h2(ii).Position(4) = originalPosition{ii}(4)*1.5;
  h2(ii).Position(2) = originalPosition{ii}(2)-0.05;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.02;
  
  h2(ii).FontSize = 12;
  colormap(h2(11),cmap)
end

for ii = [6:10 16:20]
  h2(ii).Position(2) = h2(ii).Position(2)+0.026;
end
  
for ii = 11:20
  h2(ii).Position(2) = h2(ii).Position(2)-0.05;
end

%% Figure 4: For paper, 1 sc, horizontal reversed time positioning, only perpendicular plane
scList = [1];
nsc = numel(scList);  

% Initialize particle distribution plot
nRows = 1;
nCols = 10;
isub = 1;
clear h2;

for ii = 1:(nRows*nCols);
  h2(isub) = subplot(nRows,nCols,ii); 
  h2(isub).Position(3) = h2(isub).Position(3)*1.2;
  isub = isub + 1; 
end


% Decide times
indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
indTimes = {[1214 1212 1210 1208 1206 1204],[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
indTimes = {[1215:-1:1100],[1209:-1:1100],[1218:-1:1100]};
isub = 1;
for iicc = 1:numel(scList)
  ic = scList(iicc);
  c_eval('dist = ePDist?;',ic)
  times = dist.time;
  indTime = indTimes{ic};
  % Projection coordinate system
  csys = 1;
  switch csys
    case 1 % z: B, x: N, y: zxx    
      c_eval('z = gseB?.resample(ePDist?);',ic); 
      z = -z/z.abs;
      tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
      tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

      x = cross(z,cross(tsN,z));    
      y1 = cross(z,x);
      y2 = cross(z,cross(tsL,z));    
      y = y2;
      vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
    case 2 % ExB, 3rd, B
      x = hatExB0;
      y = cross(z,x);
      z = hatB0;
      vlabels = {'B','E\times B','B\times(E\times B)'};
    case 3
      x = [1 0 0];
      y = [0 1 0];
      z = [0 0 1];
      vlabels = {'X','Y','Z'};
    case 4
      x = L;
      y = M;
      z = N;
      vlabels = {'L','M','N'};
  end
  X = x;
  Y = y;
  Z = z;

  c_eval('dist = ePDist?;',ic)
  c_eval('scpot = scPot?;',ic)
  vlim = 15*1e3;
  elevlim = 15;
  strCMap = 'jet';
  projclim = [0 4.8];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  haveYLabel = 0;
  for ii = 1:nCols % plot mms1, plane: NxB, N 
    x = X(indTime(ii)).data;
    y = Y(indTime(ii)).data;
    z = Z(indTime(ii)).data;

    time = times(indTime(ii));
    timeUTC = time.utc;      

    % Get mean vectors
    c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
    hatVe0 = double(irf_norm(Ve0));    
    c_eval('B0 = gseB?.resample(time).data;',ic); 
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = gseE?.resample(time).data;',ic); 
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);
   % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};


    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);    
    %hca.Title.String = timeUTC(15:23);
    hca.Title.String = '';
    hca.Title.String = timeUTC(15:23);
    ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    ht.VerticalAlignment = 'top';
    ht.HorizontalAlignment = 'left';
    ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    %hca.XAxisLocation = 'top';


    if ~haveYLabel          
      hca.YLabel.String = 'N_{\perp}\times B';    
      haveYLabel = 1;
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end
    %hca.XTickLabel = ''; 
    %hca.XLabel.String = ''; 
    hca.XLabel.String = 'N_{\perp}';    
  end   

end
for ii = 1:nRows*nCols
  originalPosition{ii} = h2(ii).Position;
end

hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.06;

%% Figure 4: For paper, 1 sc, horizontal reversed time positioning (or not)
ic = 1;
zdists = {[7 6.5 5.5 4.5 3.5 3 1.5 0.1 0.1]};
indTimes = {[1213 1212 1211 1210 1209 1208 1207 1206 1205]};
zdists = {zdists{ic}(1:2:end)}; indTimes = {indTimes{ic}(1:2:end)};
% Initialize particle distribution plot
nRows = 4;
nCols = 5;
isub = 0;
clear h2;

for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub + 1;         
    h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
  end
end

% Decide times
%indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
%indTimes = {[1214 1212 1210 1208 1206 1204],[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};

isub = 1;

c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?;',ic)
times = dist.time;
indTime = indTimes{ic};

% Projection coordinate system
csys = 1;
switch csys
  case 1 % z: B, x: N, y: zxx    
    c_eval('z = gseB?.resample(ePDist?);',ic); 
    z = -z/z.abs;
    tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
    tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

    x = cross(z,cross(tsN,z));    
    y1 = cross(z,x);
    y2 = cross(z,cross(tsL,z));    
    y = y1;
    vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  case 2 % ExB, 3rd, B
    x = hatExB0;
    y = cross(z,x);
    z = hatB0;
    vlabels = {'B','E\times B','B\times(E\times B)'};
  case 3
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 4
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end
X = x;
Y = y;
Z = z;

vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

haveYLabel = 0;
for ii = 1:nCols % loop over times
  % plot mms1, plane: NxB, N 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;

  time = times(indTime(ii));
  timeUTC = time.utc;      

  % Get mean vectors
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};
   
  if nRows > 0 % Perpendicular plane    
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'N\times B';                    
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'N_{\perp}';   
  end
  if nRows > 1 % Perp plane, but with markings and data of partial moments
    c_eval('neCore = ne?Core.resample(time).data;',ic)
    c_eval('neCrescent = ne?Crescent.resample(time).data;',ic)
    c_eval('veCore = ve?Core.resample(time).data;',ic)
    c_eval('veCrescent = ve?Crescent.resample(time).data;',ic)
    vectors = {hatExB0,'ExB'; E0,'E'}; % veCrescent,['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'];

    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vectors',vectors);        
    hca.Title.String = '';

    % Mark velocities and different distributions on plot
    hold(hca,'on')

    c_eval('oce = fce?.resample(time).data*2*pi;',ic)
    c_eval('ExB = gseVExB?.resample(time).data;',ic)
    c_eval('E = gseE?perp.resample(time).abs;',ic)
    c_eval('B = gseB?.resample(time).abs;',ic)     

    c_eval('v0 = minV?.resample(time).data;',ic) % km/s
    % Core 
    angles = 0:360;
    %v0 = 4.5; % 10^3 km/s
    x0 = v0*cosd(angles);
    y0 = v0*sind(angles);
    h_vcirc = plot(hca,x0*1e-3,y0*1e-3,'k');

     % Outside accelerated edge
    angles = 0:360;
    v1 = 8.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    plot(hca,x1,y1,'k');
    
    angles = 0:360;
    v1 = 7.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    plot(hca,x1,y1,'k');

    % 'vperp' lim
    zdist = 6; % km
    zdist = zdists{ic}(ii);
    vz = -10e3:1e3:10e3;
    vy = vz.^2/oce/zdist-norm(ExB)-0.25*oce*zdist;
    ang =  -15; -90+acosd(x*irf_norm(Ve0)');
    newvz = vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
    newvy = -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
    h_vcrescent = plot(hca,newvz,newvy,'k');

    % Add some info
    if 0 % separate spots
      ht = text(hca.XLim(2),hca.YLim(2),['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],'color',[1 1 1]);
      ht.VerticalAlignment = 'top';
      ht.HorizontalAlignment = 'left';
      ht.FontSize = 12; 

      htd = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'top';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 12;          

      ht = text(hca.XLim(1),hca.YLim(1),['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'],'color',[1 1 1]);
      ht.VerticalAlignment = 'bottom';
      ht.HorizontalAlignment = 'right';
      ht.FontSize = 12;        
    elseif 1 % all gathered

      if 0
        ht = text(hca.XLim(1),hca.YLim(2),{['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s']},'color',[1 1 1]);
        ht.VerticalAlignment = 'top';
        ht.HorizontalAlignment = 'right';
        ht.FontSize = 9; 
      end
      %ht = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'right';
      %ht.FontSize = 12;   

      htd = text(hca.XLim(1),hca.YLim(1),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'bottom';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 10;       
      htd.FontWeight = 'bold'; 
    end

    hold(hca,'off')


    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'N\times B';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'N_{\perp}';   
  end
  if nRows > 2 % B plane
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[x;z;-y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
      hca.Title.String = '';
      %hca.Title.String = timeUTC(15:23);
      %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'left';
      %ht.FontSize = 12;
      colormap(hca,strCMap)
      hca.XDir = 'reverse';
      if ii == 1;          
        hca.YLabel.String = 'B';          
      else
        hca.YLabel.String = '';
        hca.YTickLabel = ''; 
      end    
      hca.XLabel.String = 'N_{\perp}';   
    end
  if nRows > 3 % B plane
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    %hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      hca.YLabel.String = 'B';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    %hca.XLabel.String = 'N_{\perp}\times B';   
    hca.XLabel.String = 'N\times B';   
  end
  end   


for ii = 1:(nRows*nCols)
  originalPosition{ii} = h2(ii).Position;  
end

for ii = 1:(nRows*nCols)
  h2(ii).Position = originalPosition{ii};  
end
for ii = 1:(nRows*nCols)
  h2(ii).Position(3) = originalPosition{ii}(3)+0.1;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.05;
end

for ii = 1:4:(nRows*nCols)%[1  5  9  13  17]
  h2(ii).Position(2) = originalPosition{ii}(2);
  h2(ii).XTick = [];
  h2(ii).XLabel.String = ''; 
end
for ii = 2:4:(nRows*nCols)
  h2(ii).Position(2) = originalPosition{ii}(2)+0.05;
  h2(ii).XTick = [];
  h2(ii).XLabel.String = '';
end
for ii = 3:4:(nRows*nCols)
  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
end

for ii = 4:4:(nRows*nCols) 
  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
end

%for ii = [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19]
%  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
%end
%
hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.01;

delete(h_vcrescent)
delete(htd)

legRows = {'a','b','c','d','e','f'};
legCols = {'1)','2)','3)','4)','5)','6)'};
isub = 0;
for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub +1; % jj+(ii-1)*nCols;
    irf_legend(h2(isub),[legRows{ii} legCols{jj}],[0.02 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
end

%% Figure 4: For paper, 1 sc, first time to the left, 'marked' plot moved to the bottom 
% run first mms_2015Oct16.partial_moments
ic = 1;
zdists = {sort([7 6.5 5.5 4.5 3.5 3 1.5 0.1 0.1])};
indTimes = {sort([1213 1212 1211 1210 1209 1208 1207 1206 1205])-5*3*0};-103

step = 2;
zdists = {zdists{ic}(1:step:end)}; indTimes = {indTimes{ic}(1:step:end)};

% Initialize particle distribution plot
nRows = 4;
nCols = 5;
isub = 0;
clear h2;

for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub + 1;         
    h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
  end
end

% Decide times
%indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
%indTimes = {[1214 1212 1210 1208 1206 1204],[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};

isub = 1;

c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
times = dist.time;
indTime = indTimes{ic};

% Projection coordinate system
csys = 1;
switch csys
  case 1 % z: B, x: N, y: zxx    
    c_eval('z = gseB?.resample(ePDist?);',ic); 
    z = -z/z.abs;
    tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
    tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

    x = cross(z,cross(tsN,z));    
    y1 = cross(z,x);
    y2 = cross(z,cross(tsL,z));    
    y = y1;
    vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  case 2 % ExB, 3rd, B
    x = hatExB0;
    y = cross(z,x);
    z = hatB0;
    vlabels = {'B','E\times B','B\times(E\times B)'};
  case 3
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 4
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end
X = x;
Y = y;
Z = z;

vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

haveYLabel = 0;
for ii = 1:nCols % loop over times
  % plot mms1, plane: NxB, N 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;

  time = times(indTime(ii));
  timeUTC = time.utc;      

  % Get mean vectors
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};
   
  if nRows > 0 % Perpendicular plane    
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'v_{N\times B} (10^3 km/s)';                    
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
  end
  if nRows > 1 % B plane
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[x;z;-y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
      hca.Title.String = '';
      %hca.Title.String = timeUTC(15:23);
      %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'left';
      %ht.FontSize = 12;
      colormap(hca,strCMap)
      hca.XDir = 'reverse';
      if ii == 1;          
        hca.YLabel.String = 'v_{B} (10^3 km/s)';          
      else
        hca.YLabel.String = '';
        hca.YTickLabel = ''; 
      end    
      hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
    end
  if nRows > 2 % B plane
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    %hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      hca.YLabel.String = 'v_{B} (10^3 km/s)';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    %hca.XLabel.String = 'N_{\perp}\times B';   
    hca.XLabel.String = 'v_{N\times B} (10^3 km/s)';   
  end
  if nRows > 3 % Perp plane, but with markings and data of partial moments
    c_eval('neCore = ne?Core.resample(time).data;',ic)
    c_eval('neCrescent = ne?Crescent.resample(time).data;',ic)
    c_eval('veCore = ve?Core.resample(time).data;',ic)
    c_eval('veCrescent = ve?Crescent.resample(time).data;',ic)
    vectors = {hatExB0,'ExB'; E0,'E'}; % veCrescent,['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'];

    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vectors',vectors);        
    hca.Title.String = '';

    % Mark velocities and different distributions on plot
    hold(hca,'on')

    c_eval('oce = fce?.resample(time).data*2*pi;',ic)
    c_eval('ExB = gseVExB?.resample(time).data;',ic)
    c_eval('E = gseE?perp.resample(time).abs;',ic)
    c_eval('B = gseB?.resample(time).abs;',ic)     

    c_eval('v0 = minV?.resample(time).data;',ic) % km/s
    % Core 
    angles = 0:360;
    %v0 = 5.5*1e3; % 10^3 km/s
    x0 = v0*cosd(angles);
    y0 = v0*sind(angles);
    h_vcirc = plot(hca,x0*1e-3,y0*1e-3,'k');

     % Outside accelerated edge
    angles = 0:360;
    v1 = 8.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');
    
    angles = 0:360;
    v1 = 7.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');

    % 'vperp' lim
    zdist = 6; % km
    zdist = zdists{ic}(ii);
    vz = -10e3:1e3:10e3;
    vy = vz.^2/oce/zdist-norm(ExB)-0.25*oce*zdist;
    ang =  -15; -90+acosd(x*irf_norm(Ve0)');
    newvz = vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
    newvy = -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
    h_vcrescent = plot(hca,newvz,newvy,'k');
    if ii == 1
      delete(h_vcrescent)
    end

    % Add some info
    if 0 % separate spots
      ht = text(hca.XLim(2),hca.YLim(2),['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],'color',[1 1 1]);
      ht.VerticalAlignment = 'top';
      ht.HorizontalAlignment = 'left';
      ht.FontSize = 12; 

      htd = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'top';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 12;          

      ht = text(hca.XLim(1),hca.YLim(1),['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'],'color',[1 1 1]);
      ht.VerticalAlignment = 'bottom';
      ht.HorizontalAlignment = 'right';
      ht.FontSize = 12;        
    elseif 1 % all gathered

      if 0
        ht = text(hca.XLim(1),hca.YLim(2),{['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s']},'color',[1 1 1]);
        ht.VerticalAlignment = 'top';
        ht.HorizontalAlignment = 'right';
        ht.FontSize = 9; 
      end
      %ht = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'right';
      %ht.FontSize = 12;   

      if 0
        htd = text(hca.XLim(1),hca.YLim(1),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
        htd.VerticalAlignment = 'bottom';
        htd.HorizontalAlignment = 'right';
        htd.FontSize = 10;       
        htd.FontWeight = 'bold'; 
      else % put it int he title
        hca.Title.String = ['d = ' num2str(zdist,'%g') ' km'];
        %htd.FontSize = 10;       
        hca.Title.FontWeight = 'light'; 
        if ii ==1
          hca.Title.String = '';
        end
      end
    end

    hold(hca,'off')


    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'v_{N\times B} (10^3 km/s)';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
  end
  end   


for ii = 1:(nRows*nCols)
  originalPosition{ii} = h2(ii).Position;  
end
%%
for ii = 1:(nRows*nCols)
  h2(ii).Position = originalPosition{ii};  
  h2(ii).FontSize = 12;
  h2(ii).XGrid = 'off';
  h2(ii).YGrid = 'off';
end

for ii = 1:(nRows*nCols)
  h2(ii).Position(3) = originalPosition{ii}(3)+0.1;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.05;
end

for ii = 1:nRows:(nRows*nCols) % 1st row
  h2(ii).Position(2) = originalPosition{ii}(2);
  h2(ii).XTick = [];
  h2(ii).XLabel.String = ''; 
end
for ii = 2:nRows:(nRows*nCols) % 2nd row
  h2(ii).Position(2) = originalPosition{ii}(2)+0.05;
  %h2(ii).XTick = [];
  %h2(ii).XLabel.String = '';
end
for ii = 3:nRows:(nRows*nCols) % 3rd row
  h2(ii).Position(2) = originalPosition{ii}(2)+0.04;
end

for ii = 4:4:(nRows*nCols) % 4th row
  h2(ii).Position(2) = originalPosition{ii}(2)+0.00;
end

%for ii = [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19]
%  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
%end
%
hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.01;

%delete(h_vcrescent)
%delete(htd)

legRows = {'a','b','c','d','e','f','g','h','i','j'};
legCols = {'1)','2)','3)','4)','5)','6)'};
tsshift = 3;
isub = 0;
for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub +1; % jj+(ii-1)*nCols;
    hhh = irf_legend(h2(isub),[legRows{ii+3} legCols{jj}],[0.02 0.98],'color',[0 0 0],'fontweight','bold','fontsize',12);
    hhh.BackgroundColor = [1 1 1];
    hhh.EdgeColor = [0 0 0];
  end
end
noXlabel = [1:8 13:20]; for il = noXlabel; h2(il).XLabel.String = ''; end


%% Figure 4: For paper, 1 sc, including PA 0 90 180 plot, horizontal reversed time positioning (or not)
ic = 1;
zdists = {[7 6.5 5.5 4.5 3.5 3 1.5 0.1 0.1]};
indTimes = {[1213 1212 1211 1210 1209 1208 1207 1206 1205]};
zdists = {zdists{ic}(1:2:end)}; indTimes = {indTimes{ic}(1:2:end)};
% Initialize particle distribution plot
nRows = 5;
nCols = 5;
isub = 0;
clear h2;

for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub + 1;         
    h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
  end
end

% Decide times
%indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
%indTimes = {[1214 1212 1210 1208 1206 1204],[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};

isub = 1;

c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?;',ic)
times = dist.time;
indTime = indTimes{ic};

% Projection coordinate system
csys = 1;
switch csys
  case 1 % z: B, x: N, y: zxx    
    c_eval('z = gseB?.resample(ePDist?);',ic); 
    z = -z/z.abs;
    tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
    tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

    x = cross(z,cross(tsN,z));    
    y1 = cross(z,x);
    y2 = cross(z,cross(tsL,z));    
    y = y1;
    vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  case 2 % ExB, 3rd, B
    x = hatExB0;
    y = cross(z,x);
    z = hatB0;
    vlabels = {'B','E\times B','B\times(E\times B)'};
  case 3
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 4
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end
X = x;
Y = y;
Z = z;

vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

haveYLabel = 0;
for ii = 1:nCols % loop over times
  % plot mms1, plane: NxB, N 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;

  time = times(indTime(ii));
  timeUTC = time.utc;      

  % Get mean vectors
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};
   
  if nRows > 0 % Perpendicular plane    
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'N\times B';                    
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'N_{\perp}';   
  end
  if nRows > 1 % Perp plane, but with markings and data of partial moments
    c_eval('neCore = ne?Core.resample(time).data;',ic)
    c_eval('neCrescent = ne?Crescent.resample(time).data;',ic)
    c_eval('veCore = ve?Core.resample(time).data;',ic)
    c_eval('veCrescent = ve?Crescent.resample(time).data;',ic)
    vectors = {hatExB0,'ExB'; E0,'E'}; % veCrescent,['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'];

    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vectors',vectors);        
    hca.Title.String = '';

    % Mark velocities and different distributions on plot
    hold(hca,'on')

    c_eval('oce = fce?.resample(time).data*2*pi;',ic)
    c_eval('ExB = gseVExB?.resample(time).data;',ic)
    c_eval('E = gseE?perp.resample(time).abs;',ic)
    c_eval('B = gseB?.resample(time).abs;',ic)     

    c_eval('v0 = minV?.resample(time).data;',ic) % km/s
    % Core 
    angles = 0:360;
    %v0 = 4.5; % 10^3 km/s
    x0 = v0*cosd(angles);
    y0 = v0*sind(angles);
    h_vcirc = plot(hca,x0*1e-3,y0*1e-3,'k');

     % Outside accelerated edge
    angles = 0:360;
    v1 = 8.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');
    
    angles = 0:360;
    v1 = 7.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');

    % 'vperp' lim
    zdist = 6; % km
    zdist = zdists{ic}(ii);
    vz = -10e3:1e3:10e3;
    vy = vz.^2/oce/zdist-norm(ExB)-0.25*oce*zdist;
    ang =  -15; -90+acosd(x*irf_norm(Ve0)');
    newvz = vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
    newvy = -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
    h_vcrescent = plot(hca,newvz,newvy,'k');

    % Add some info
    if 0 % separate spots
      ht = text(hca.XLim(2),hca.YLim(2),['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],'color',[1 1 1]);
      ht.VerticalAlignment = 'top';
      ht.HorizontalAlignment = 'left';
      ht.FontSize = 12; 

      htd = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'top';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 12;          

      ht = text(hca.XLim(1),hca.YLim(1),['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'],'color',[1 1 1]);
      ht.VerticalAlignment = 'bottom';
      ht.HorizontalAlignment = 'right';
      ht.FontSize = 12;        
    elseif 1 % all gathered

      if 0
        ht = text(hca.XLim(1),hca.YLim(2),{['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s']},'color',[1 1 1]);
        ht.VerticalAlignment = 'top';
        ht.HorizontalAlignment = 'right';
        ht.FontSize = 9; 
      end
      %ht = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'right';
      %ht.FontSize = 12;   

      htd = text(hca.XLim(1),hca.YLim(1),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'bottom';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 10;       
      htd.FontWeight = 'bold'; 
    end

    hold(hca,'off')


    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'N\times B';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'N_{\perp}';   
  end
  if nRows > 2 % B plane
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[x;z;-y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
      hca.Title.String = '';
      %hca.Title.String = timeUTC(15:23);
      %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'left';
      %ht.FontSize = 12;
      colormap(hca,strCMap)
      hca.XDir = 'reverse';
      if ii == 1;          
        hca.YLabel.String = 'B';          
      else
        hca.YLabel.String = '';
        hca.YTickLabel = ''; 
      end    
      hca.XLabel.String = 'N_{\perp}';   
    end
  if nRows > 3 % B plane
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    %hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      hca.YLabel.String = 'B';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    %hca.XLabel.String = 'N_{\perp}\times B';   
    hca.XLabel.String = 'N\times B';   
  end
  if nRows > 4 % Pitchangle plot, 0 90 180    
    unitscale = 1e30; % cm-6->km-6
    hca = h2(isub); isub = isub + 1;     
    plot(hca,ePitch1.depend{1}(indTime(ii),:),squeeze(ePitch1(indTime(ii)).data(:,:,1))*unitscale,...
             ePitch1.depend{1}(indTime(ii),:),squeeze(ePitch1(indTime(ii)).data(:,:,7))*unitscale,...
             ePitch1.depend{1}(indTime(ii),:),squeeze(ePitch1(indTime(ii)).data(:,:,13)*unitscale))
    hca.YScale = 'log';
    hca.XScale = 'log';
    hca.XLim = [1e1 2e2];
    hca.YLim = [5e-28 1e-25]*unitscale;
    hca.YTick = [1e-28 1e-27 1e-26 1e-25]*unitscale;
    axis(hca,'square')
    if ii == 1 
      hca.YLabel.String = 'f_e [s^3 km^-6]';
    else
      hca.YTickLabel = {};
    end
    %hleg = legend(hca,'\theta = 0','\theta = 90','\theta = 180');
    %hleg.box = 'off';
  end
end   

for ii = 1:(nRows*nCols)
  originalPosition{ii} = h2(ii).Position;  
end
%%

for ii = 1:(nRows*nCols)
  h2(ii).Position = originalPosition{ii};  
end
for ii = 1:(nRows*nCols)
  h2(ii).Position(3) = originalPosition{ii}(3)+0.1;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.05;
end

for ii = 1:nRows:(nRows*nCols)%[1  5  9  13  17]
  h2(ii).Position(2) = originalPosition{ii}(2);
  h2(ii).XTick = [];
  h2(ii).XLabel.String = ''; 
end
for ii = 2:nRows:(nRows*nCols)
  h2(ii).Position(2) = originalPosition{ii}(2)+0.05;
  h2(ii).XTick = [];
  h2(ii).XLabel.String = '';
end
for ii = 3:nRows:(nRows*nCols)
  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
end

for ii = 4:nRows:(nRows*nCols) 
  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
end


for ii = 5:nRows:(nRows*nCols) 
  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.08;
end
%for ii = [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19]
%  h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
%end
%
hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.01;

delete(h_vcrescent)
delete(htd)

legRows = {'a','b','c','d','e','f','g','h','i'};
legCols = {'1)','2)','3)','4)','5)','6)'};
nTS = 3;
isub = 0;
for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub +1; % jj+(ii-1)*nCols;
    irf_legend(h2(isub),[legRows{ii}+nTS legCols{jj}],[0.02 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
end
%% Only one time 0 90 180 plot, where the waves are observed
unitscale = 1e30; % cm-6->km-6
its = [1213 1212 1211 1210 1209 1208 1207 1206 1205];
hca = subplot(1,1,1); isub = isub + 1;    
ii = 4;
plot(hca,ePitch1.depend{1}(its(ii),:),squeeze(ePitch1(its(ii)).data(:,:,1))*unitscale,...
         ePitch1.depend{1}(its(ii),:),squeeze(ePitch1(its(ii)).data(:,:,7))*unitscale,...
         ePitch1.depend{1}(its(ii),:),squeeze(ePitch1(its(ii)).data(:,:,13)*unitscale))
hca.YScale = 'log';
hca.XScale = 'log';
hca.XLim = [1e1 5e2];
hca.YLabel.String = 'f_e [s^3km^{-6}]';
hca.XLabel.String = 'E [eV]';
hca.YLim = [1e-29 1e-25]*unitscale;
hca.YTick = [1e-28 1e-27 1e-26 1e-25]*unitscale;
axis(hca,'square')
time = ePitch1(its(ii)).time;
timeUTC = time.utc; 
hca.Title.String = timeUTC(12:23)
hleg = legend(hca,'\theta = 0','\theta = 90','\theta = 180');
hleg.EdgeColor = [1 1 1];

%% Parallel pitch angle plot, to relate to Epar waves
hca = subplot(1,1,1);
hold(hca,'on')
unitscale = 1e30;
its = [1213 1212 1211 1210 1209 1208 1207 1206 1205];
for it = its(1:8);[1213 1212 1211 1210 1209 1208 1207 1206 1205];  
  plot(hca,ePitch1.depend{1}(it,:),squeeze(ePitch1(it).data(:,:,7))*unitscale)
end
hold(hca,'off')
hca.YScale = 'log';
hca.XScale = 'log';
%hca.XLim = [1e1 2e2];
%hca.YLim = [5e-28 1e-25]*unitscale;
hca.YTick = [1e-28 1e-27 1e-26 1e-25]*unitscale;
%axis(hca,'square')
    if ii == 1 
      hca.YLabel.String = 'f_e [s^3 km^-6]';
    else
      hca.YTickLabel = {};
    end
%% Support plot to the electron distributions
tintZoom = irf.tint('2015-10-16T10:33:30.00Z/2015-10-16T10:33:30.60Z');
ic = 1;

h = irf_plot(5);

if 0 % LMN, 2 panels
  hca = irf_panel('ve bulk');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('lines_bulk = irf_plot(hca,{mvaVe?.x,mvaVe?.y,mvaVe?.z},''comp'');',ic)
  hca.YLabel.String = {'v_e','km/s'};
  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.95])
  c_eval('lines_core = irf_plot(hca,{mvaVe?Core.x,mvaVe?Core.y,mvaVe?Core.z},''comp'',''-.'');',ic)%,''.''

  h_lines = findall(hca,'type','line'); 
  h_leg = legend([h_lines([1 4]+2)],{'core','bulk'});
  h_leg.Box = 'off';

  hca = irf_panel('ve crescent');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('lines_crescent = irf_plot(hca,{mvaVe?Crescent.x,mvaVe?Crescent.y,mvaVe?Crescent.z},''comp'',''-.'')',ic)%,''.''
  set(hca,'ColorOrder',mms_colors('xyza'))
  hca.YLabel.String = {'v_e','km/s'};
  hca.YLabel.Interpreter = 'tex';
else % L, 1 panel
  hca = irf_panel('ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('lines_bulk = irf_plot(hca,{mvaVe?.x,mvaVe?Core.x,mvaVe?Crescent.x},''comp'');',ic)
  hca.YLabel.String = {'v_{e,L}','km/s'};
  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'Total','Core','Crescent'},[0.98 0.5])
end
%c_eval('lines_crescent = irf_plot(hca,{mvaVe?Core.x,mvaVe?Core.y,mvaVe?Core.z},''comp'','':'')',ic)%,''.''

hca = irf_panel('ne');
set(hca,'ColorOrder',mms_colors('1234b'))
c_eval('irf_plot(hca,{ne?,ne?Core,ne?Crescent},''comp'');',ic)
%c_eval('irf_plot(hca,{ne?,ne?Core,ne?Crescent,ne?BG,ne?Crescent-ne?BG},''comp'');',ic)
hca.YLabel.String = {'n_e','cm^{-3}'};
set(hca,'ColorOrder',mms_colors('1234b'))
%irf_legend(hca,{'total','Core','Crescent'},[0.98 0.95])

hca = irf_panel('Q');
set(hca,'ColorOrder',mms_colors('1234'))
c_eval('irf_plot(hca,{sqrt(Q?),sqrt(Q?Core),sqrt(Q?Crescent)},''comp'');',ic)
%hca.YLabel.String = {'$$\sqrt{Q}$$',''};
%hca.Interpreter = 'latex';
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'total','Core','Crescent'},[0.98 0.95])
ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');

if 0
  hca = irf_panel('Ao');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{Ao?,Ao?Core,Ao?Crescent},''comp'');',ic)
  %hca.YLabel.String = {'$$\sqrt{Q}$$',''};
  %hca.Interpreter = 'latex';
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'total','Core','Crescent'},[0.98 0.95])
  ylabel(hca,'Ao','Interpreter','latex');
end
if 1
  hca = irf_panel('Tratio');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{T?ratio,T?ratioCore,T?ratioCrescent},''comp'');',ic)
  %hca.YLabel.String = {'$$\sqrt{Q}$$',''};
  %hca.Interpreter = 'latex';
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'total','Core','Crescent'},[0.98 0.95])
  ylabel(hca,'T_{||}/T_{\perp}','Interpreter','tex');
end
if 1 % Epar
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1par.tlim(tint),mvaE2par.tlim(tint),mvaE3par.tlim(tint),mvaE4par.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:numel(h)
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
end
for ii=1:numel(indTime)
  %c_eval('irf_pl_mark(h,ePDist?.time(indTime(ii)).epochUnix,[0.5 0.5 0.5]);',ic)
  c_eval('irf_pl_mark(h,ePDist?.time(indTime(ii)).epochUnix,[0.5 0.5 0.9]);',ic)
end

%% Figure 4: For paper, 1 sc
ic = 1;
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');

% Initialize timeseries plot
h1 = irf_plot(1);
h1.Position(2) = 0.85; h1.Position(4) = 0.08;

hca = irf_panel('Ve');
c_eval('irf_plot(hca,{mvaVe?par,mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)
irf_zoom(h1(1),'x',tintZoom)
irf_zoom(h1(1),'y')
irf_plot_axis_align

% Initialize particle distribution plot
nrows = 3;
ncols = 4;
isub = 1;
for ii = ncols+1:(nrows+1)*ncols; h2(isub) = subplot(nrows+1,ncols,ii); isub = isub + 1; end
xLeft = h2(1).Position(1); xRight = h2(4).Position(1)+h2(4).Position(3);
xGap = 0.05;
xWidth = (xRight-xLeft)/ncols - (ncols/2-1)*xGap;

yTop = h2(1).Position(2)+h2(1).Position(3); yBot = h2(end).Position(2);
yGap = 0.05;
yWidth = (yTop-yBot)/nrows - (nrows/2-1)*yGap;
xWidth = yWidth;

ii = 1;
for irow = 1:nrows
  for icol = 1:ncols          
    xgap = ceil(icol/2-1)*xGap;
    x0 = xLeft + xgap + (icol-1)*xWidth;
    ygap = (irow-1)*yGap;
    y0 = yTop - yWidth - ygap - (irow-1)*yWidth;
        
    pos = [x0 y0 xWidth yWidth];
    h2(ii).Position = pos;
    ii = ii + 1;    
  end
end

% Decide times
c_eval('dist = ePDist?;',ic)
times = dist.time;
indTime = 1210+[0:2:11]-5;

% Projection coordinate system
csys = 1;
switch csys
  case 1 % z: B, x: N, y: zxx    
    c_eval('z = gseB?.resample(ePDist?);',ic); 
    z = -z/z.abs;
    tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
    tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));
    x = cross(z,cross(tsN,z));    
    y = cross(z,x);       
    vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  case 2
    x = hatB0;
    y = hatExB0;
    z = cross(x,y);
    vlabels = {'B','E\times B','B\times(E\times B)'};
  case 3
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 4
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end
X = x;
Y = y;
Z = z;

c_eval('dist = ePDist?;',ic)
vlim = 15*1e3;
elevlim = 25;
strCMap = 'jet';
%energies =  [30 220];
projclim = [0 4.5];
palim = [1e-3 1e6];
skymapEnergy = [65 278];


isub = 1;
for ii = 1:numel(indTime) 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;
  
  time = times(indTime(ii));
  
  hmark = irf_pl_mark(h1(1),[time.epochUnix],'yellow');
  
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    

  % Get mean magnetic field direction
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);

 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

  
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  hca.XDir = 'reverse';
  
  hca.XLabel.String = '';
  hca.YLabel.String = '';
  hca.XLabel.String = 'N_{\perp}';
  hca.YLabel.String = 'N_{\perp}\times B';
  
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  mms.plot_projection(hca,dist,'tint',time,'xyz',[-y;-z;-x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  hca.XDir = 'reverse';
  
  hca.XLabel.String = 'N_{\perp}';
  hca.YLabel.String = 'B';
end 
hcf = gcf;

cb = hcf.Children(24);
delete(hcf.Children(2:2:22))

cb.Position = [h2(ncols).Position(1)+xWidth+0.03 h2(end).Position(2)+0.01 0.03 yTop-yBot];

% Put labels on just some axes
hca.XLabel.String = vlabels{1};
hca.YLabel.String = vlabels{2};

%% Figure 4: Particle distributions
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h1 = irf_plot(8);
isub = 1;
for ii = 4:15; h2(isub) = subplot(5,3,ii); isub = isub + 1;end

for ic = 1
  hca = irf_panel(irf_ssub('Ve mms ?',ic));
  set(hca,'ColorOrder',mms_colors('axyz'))
  c_eval('irf_plot(hca,{mvaVe?par,mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('axyz'))
  irf_legend(hca,{'||','\perp_L','\perp_M','\perp_N'},[0.02 0.9],'fontsize',12);
end

ic = 1;
c_eval('dist = ePDist?;',ic)
times = dist.time;
indTime = [1209 1210 1211 1212]; % +10;

vlim = 15*1e3;
elevlim = 25;
strCMap = 'jet';
%energies =  [30 220];
projclim = [0 4.5];
palim = [1e-3 1e6];
skymapEnergy = [65 278];

isub = 1;
for ii = 1:numel(indTime)
  time = times(indTime(ii));
  hmark = irf_pl_mark(h1(1),time.epochUnix','green');
  
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    

  % Get mean magnetic field direction
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);

  vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

  % Projection coordinate system
  if 1
    x = hatB0;
    y = hatExB0;
    z = cross(x,y);
  else
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
  end
%  isub = 1;

  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)

end


irf_zoom(h1(1),'x',tintZoom)
irf_zoom(h1(1),'y')
irf_plot_axis_align

%% Figure 4: Particle distributions only perpendicular plane
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h1 = irf_plot(8);
isub = 1;
for ii = 4:15; h2(isub) = subplot(5,3,ii); isub = isub + 1;end
sclist = 1;
for ic = sclist
  hca = irf_panel(irf_ssub('Ve mms ?',ic));
  set(hca,'ColorOrder',mms_colors('axyz'))
  c_eval('irf_plot(hca,{mvaVe?par,mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('axyz'))
  irf_legend(hca,{'||','\perp_L','\perp_M','\perp_N'},[0.02 0.9],'fontsize',12);
end

ic = sclist;
c_eval('dist = ePDist?;',ic)
c_eval('scpot = scPot?;',ic)
times = dist.time;
indTime = [1209 1210 1211 1212];
indTime = 1207+[0:11]-5+10;

vlim = 15*1e3;
elevlim = 25;
strCMap = 'jet';
%energies =  [30 220];
projclim = [0 4.5];
palim = [1e-3 1e6];
skymapEnergy = [65 278];
hmark = irf_pl_mark(h1(1),[times(indTime(1)).epochUnix times(indTime(end)).epochUnix],'yellow');

isub = 1;
for ii = 1:numel(indTime)
  time = times(indTime(ii));
  
  
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    

  % Get mean magnetic field direction
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);

  vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

  % Projection coordinate system
  if 1
    z = -hatB0;
    x = irf_norm(cross(hatB0,cross(N,hatB0)));
    y = cross(z,x);   
    vlabels = {'B\times(N\times B)','3d','B'}
  elseif 1
    x = hatB0;
    y = hatExB0;
    z = cross(x,y);
  else
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
  end
  
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'vlabel',vlabels,'flipx','scpot',scpot);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)

end


irf_zoom(h1(1),'x',tintZoom)
irf_zoom(h1(1),'y')
irf_plot_axis_align

%% Figure: Electron moments across boundary, time shifted
[h,h2] = initialize_combined_plot(8,1,2,0.3,'vertical');
dt = [0 -0.07 -0.20 -0.02]*0;

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z'); 
tintZoom = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:32.00Z'); 

ic = 1;
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni 32');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni 64');
  c_eval('irf_spectrogram(hca,ePDist?.e64.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 0 % ePDist pa 32 0 200
  hca = irf_panel('e PA deflux low');  
  eint = [0 200];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,[]).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  irf_legend(hca,['E: ' num2str(eint(1),'%g') '-' num2str(eint(2),'%g') ' eV'],[0.95 0.03])
end
if 0 % ePDist pa 32: 200 40000
  hca = irf_panel('e PA deflux high');  
  eint = [200 40000];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,[]).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  irf_legend(hca,['E: ' num2str(eint(1),'%g') '-' num2str(eint(2),'%g') ' eV'],[0.95 0.03])
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp','dt',dt);    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp','dt',dt);    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp','dt',dt);    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp','dt',dt);    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 0 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).x,mvaE2.tlim(tint).x,mvaE3.tlim(tint).x,mvaE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).y,mvaE2.tlim(tint).y,mvaE3.tlim(tint).y,mvaE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EL par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1par.tlim(tint),mvaE2par.tlim(tint),mvaE3par.tlim(tint),mvaE4par.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EL perp
  hca = irf_panel('EL perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).x,mvaE2perp.tlim(tint).x,mvaE3perp.tlim(tint).x,mvaE4perp.tlim(tint).x},'comp','dt',dt);
  hca.YLabel.String = {'E_{\perp,L}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM perp
  hca = irf_panel('EM perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).y,mvaE2perp.tlim(tint).y,mvaE3perp.tlim(tint).y,mvaE4perp.tlim(tint).y},'comp','dt',dt);
  hca.YLabel.String = {'E_{\perp,M}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EN perp
  hca = irf_panel('EN perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).z,mvaE2perp.tlim(tint).z,mvaE3perp.tlim(tint).z,mvaE4perp.tlim(tint).z},'comp','dt',dt);
  hca.YLabel.String = {'E_{\perp,N}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0
%load('caa/cmap.mat');
%for ii = [5 6]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
%  colormap(h(ii),'jet')
%end
%legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
%for ii = 1:10
%  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
%end
end


irf_zoom(h,'x',tintZoom)
%irf_zoom(h([1:3 6 10]),'y')
irf_plot_axis_align
%

%% Plot quivers in right panels

isub = 1;      
tintQuivers = irf.tint('2015-10-16T10:33:29.30Z/2015-10-16T10:33:31.10Z');
if 0 % spacecraft configuration   
  hca = h2(isub); isub = isub + 1; 
  lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 -1 0])
  axis(hca,'square')
end
if 0
  hca = h2(isub); isub = isub +1;
  hold(hca,'on')
  irf_pl_mark(h,tintQuivers.epochUnix')
  times = mvaVe1.tlim(tintQuivers).time;
  gseV = 55*[-0.90 -0.28 -0.33]; % GSE
  lmnV = gseV*[L' M' N'];
  c_eval('posR? = repmat(mvaRR?,times.length,1)-(times-times(1))*lmnV;')
  c_eval('posV? = mvaVe?.resample(times).data;')
  c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  hold(hca,'off')

  hca.YLabel.String = 'L';
  hca.XLabel.String = 'N';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0
  hca = h2(isub); isub = isub +1;
  hold(hca,'on')
  irf_pl_mark(h,tintQuivers.epochUnix')
  times = mvaVe1.tlim(tintQuivers).time;
  gseV = 55*[-0.90 -0.28 -0.33]; % GSE
  tanV = irf_norm(cross(gseV,M));
  gseV = 55*[-0.90 -0.28 -0.33]-20*tanV;
  lmnV = gseV*[L' M' N'];
  c_eval('posR? = repmat(mvaRR?,times.length,1)-(times-times(1))*lmnV;')
  c_eval('posV? = mvaVe?.resample(times).data;')
  c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  hold(hca,'off')

  hca.YLabel.String = 'L';
  hca.XLabel.String = 'N';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % 3D
 
  irf_pl_mark(h,tintQuivers.epochUnix')
  times = mvaVe1.tlim(tintQuivers).time;
  
  gseVouteredge = 65.3*[0.37  0.32 -0.87];  
  gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
  gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
  gseMSP = 32*[-0.88 -0.42 0.24];
  
  vel_selection = 2;
  clear timesVUTC gseVdata
  switch vel_selection
    case 1 % do velocities manually    
      tanV = irf_norm(cross(gseV,M));
      gseV = 55*[-0.90 -0.28 -0.33]-30*tanV;

      gseV = repmat(gseV,times.length,1);
      gseV(1:27,:) = repmat(-50*L,27,1);   
      gseV(28:31,:) = repmat(-50*N,4,1);   
      gseV(44:60,:) = repmat(gseVouteredge,17,1);  
    case 2 % define velocities at certain times, and then interpolate to other times      
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-25*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 3 % more 'vertical'
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
  end
  
  dt = times(2)-times(1);
  
  lmnV = gseV.data*[L' M' N'];
  
  
  c_eval('posR? = repmat(mvaRR?,times.length,1)-dt*cumsum(lmnV,1);')
  c_eval('posV? = mvaVe?.resample(times).data;')
  c_eval('posJ? = mvaJ?.resample(times).data;')
  c_eval('posB? = mvaB?.resample(times).data;')
  c_eval('posE? = mvaE?.resample(times).data;')
  c_eval('posRe? = mvaEVexB?.resample(times).data;')
  
  hca = h2(isub); isub = isub +1;
  sclist = 1;
  hold(hca,'on')
  %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posV?(:,3) -posV?(:,2) posV?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
  %c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
    c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
  hold(hca,'off')

  hca.XLabel.String = 'N';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'L';
  hca.ZDir = 'normal';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  %axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  view(hca,[0 -1 0])
  %hca.YLim = [-20 20];
  %axis(hca,'equal')
  axis(hca,'equal')
  hca.Title.String = 'V_e, E (faint)';
  
  if 1 % plot electric field arrows
    hca = h2(isub); isub = isub +1;
    hold(hca,'on')
    %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
    c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
    c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
    c_eval('plot_quivers(hca,[posJ?(:,3) -posJ?(:,2) posJ?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
    %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',1:4)
    %c_eval('plot_quivers(hca,[posRe?(:,3) -posRe?(:,2) posRe?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
    hold(hca,'off')

    hca.XLabel.String = 'N';
    hca.YLabel.String = 'M';
    hca.ZLabel.String = 'L';
    hca.ZDir = 'normal';
    hca.YDir = 'normal';
    hca.XDir = 'reverse';
    %axis(hca,'square')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.ZGrid = 'on';
    view(hca,[0 -1 0])
    %hca.YLim = [-20 20];
    %axis(hca,'equal')
    %axis(hca,'equal')
    %hca.XLim = 20*[-1 1];
    hca.Title.String = 'E, J (faint)';
  end
  if 0 % Interpolate BM in LN-plane to get Hall field color surface
    posL = double([posR1(:,1); posR2(:,1); posR3(:,1); posR4(:,1)]);
    posM = double([posR1(:,2); posR2(:,2); posR3(:,2); posR4(:,2)]);
    posN = double([posR1(:,3); posR2(:,3); posR3(:,3); posR4(:,3)]);
    posB = double([posB1(:,2); posB2(:,2); posB3(:,2); posB4(:,2)]);
    dN = 2; dL = 2;
    [NN,LL] = meshgrid(min(posN):dN:max(posN),min(posL):dL:max(posL));
    fBM = griddata(posN,posL,posB,NN,LL);
    hold(hca,'on');
    mesh(hca,NN,NN*0+abs(min(posM)),LL,fBM);  
    hmcb = colorbar('peer',hca); 
    %plot3(hca,posN,posM,posL,posB,'o');
    hold(hca,'off');
  end
end

%% Quiver plot for paper
h2 = subplot(1,1,1);
isub = 1;      
tintQuivers = irf.tint('2015-10-16T10:33:29.50Z/2015-10-16T10:33:30.60Z');

if 0
  hca = h2(isub); isub = isub +1;
  hold(hca,'on')
  irf_pl_mark(h,tintQuivers.epochUnix')
  times = mvaVe1.tlim(tintQuivers).time;
  gseV = 55*[-0.90 -0.28 -0.33]; % GSE
  lmnV = gseV*[L' M' N'];
  c_eval('posR? = repmat(mvaRR?,times.length,1)-(times-times(1))*lmnV;')
  c_eval('posV? = mvaVe?.resample(times).data;')
  c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  hold(hca,'off')

  hca.YLabel.String = 'L';
  hca.XLabel.String = 'N';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0
  hca = h2(isub); isub = isub +1;
  hold(hca,'on')
  irf_pl_mark(h,tintQuivers.epochUnix')
  times = mvaVe1.tlim(tintQuivers).time;
  gseV = 55*[-0.90 -0.28 -0.33]; % GSE
  tanV = irf_norm(cross(gseV,M));
  gseV = 55*[-0.90 -0.28 -0.33]-20*tanV;
  lmnV = gseV*[L' M' N'];
  c_eval('posR? = repmat(mvaRR?,times.length,1)-(times-times(1))*lmnV;')
  c_eval('posV? = mvaVe?.resample(times).data;')
  c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  hold(hca,'off')

  hca.YLabel.String = 'L';
  hca.XLabel.String = 'N';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % 3D
   
  times = mvaVe1.tlim(tintQuivers).time;
  
  gseVouteredge = 65.3*[0.37  0.32 -0.87];  
  gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
  gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
  gseMSP = 32*[-0.88 -0.42 0.24];
  
  lmnVmsh = gseVinneredge*[L' M' N'];
  lmnVmsp = gseMSP*[L' M' N'];
  
  vel_selection = 5;
  clear timesVUTC gseVdata
  switch vel_selection
    case 1 % do velocities manually    
      tanV = irf_norm(cross(gseV,M));
      gseV = 55*[-0.90 -0.28 -0.33]-30*tanV;

      gseV = repmat(gseV,times.length,1);
      gseV(1:27,:) = repmat(-50*L,27,1);   
      gseV(28:31,:) = repmat(-50*N,4,1);   
      gseV(44:60,:) = repmat(gseVouteredge,17,1);  
    case 2 % define velocities at certain times, and then interpolate to other times      
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-25*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 3 % more 'vertical'
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 4
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 5
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge)...
                                                                                 -cross(irf_norm(gseMSP),M)*33*0.75;      
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge*0.75; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 6
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge); 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge;-40*tanVinneredge; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge;+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
  end
  
  dt = times(2)-times(1);
  
  lmnV = gseV.data*[L' M' N'];
  
  
  c_eval('posR? = repmat(mvaRR?,times.length,1)-dt*cumsum(lmnV,1);')
  c_eval('posV? = mvaVe?.resample(times).data*0.6;')
  c_eval('posJ? = mvaJ?.resample(times).data;')
  c_eval('posB? = mvaB?.resample(times).data;')
  c_eval('posE? = mvaE?.resample(times).data;')
  c_eval('posRe? = mvaEVexB?.resample(times).data;')
  
  hca = h2(isub); isub = isub +1;
  sclist = 1:4;
  hold(hca,'on')
  %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posV?(:,3) -posV?(:,2) posV?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
  %c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
  %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
  hold(hca,'off')

  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'M (km)';
  hca.ZLabel.String = 'L (km)';
  hca.ZDir = 'normal';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  %axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  view(hca,[0 -1 0])
  %hca.YLim = [-20 20];
  %axis(hca,'equal')
  axis(hca,'equal')
  
  hca.Box = 'on';
  hca.ZLim = [-10 75];
  hca.XLim = [-8 27];
  tintUTCstart = tintQuivers(1).utc;
  tintUTCstop = tintQuivers(2).utc;
  hca.Title.String = ['V_e: ' tintUTCstop(12:22) ' - ' tintUTCstop(12:22)];
  hca.Title.String = ['V_e'];
  hleg=legend(hca,'MMS 1','MMS 2','MMS 3','MMS 4','location','northeast');
  fontsize = 17;
  hca.XLabel.FontSize = fontsize;
  hca.YLabel.FontSize = fontsize;
  hca.ZLabel.FontSize = fontsize;
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  hleg.FontSize = 10;
  if 0 % plot electric field arrows
    hca = h2(isub); isub = isub +1;
    hold(hca,'on')
    %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
    c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
    c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
    c_eval('plot_quivers(hca,[posJ?(:,3) -posJ?(:,2) posJ?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
    %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',1:4)
    %c_eval('plot_quivers(hca,[posRe?(:,3) -posRe?(:,2) posRe?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
    hold(hca,'off')

    hca.XLabel.String = 'N';
    hca.YLabel.String = 'M';
    hca.ZLabel.String = 'L';
    hca.ZDir = 'normal';
    hca.YDir = 'normal';
    hca.XDir = 'reverse';
    %axis(hca,'square')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.ZGrid = 'on';
    view(hca,[0 -1 0])
    %hca.YLim = [-20 20];
    %axis(hca,'equal')
    %axis(hca,'equal')
    %hca.XLim = 20*[-1 1];
    hca.Title.String = 'E, J (faint)';
  end
  if 0 % Interpolate BM in LN-plane to get Hall field color surface
    posL = double([posR1(:,1); posR2(:,1); posR3(:,1); posR4(:,1)]);
    posM = double([posR1(:,2); posR2(:,2); posR3(:,2); posR4(:,2)]);
    posN = double([posR1(:,3); posR2(:,3); posR3(:,3); posR4(:,3)]);
    posB = double([posB1(:,2); posB2(:,2); posB3(:,2); posB4(:,2)]);
    dN = 2; dL = 2;
    [NN,LL] = meshgrid(min(posN):dN:max(posN),min(posL):dL:max(posL));
    fBM = griddata(posN,posL,posB,NN,LL);
    hold(hca,'on');
    mesh(hca,NN,NN*0+abs(min(posM)),LL,fBM);  
    hmcb = colorbar('peer',hca); 
    %plot3(hca,posN,posM,posL,posB,'o');
    hold(hca,'off');
  end
  
  
  if 1 % plot magnetosheath and magnetosphere boundary planes    
    gseVouteredge = 65.3*[0.37  0.32 -0.87];  
    gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
    gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
    gseMSP = 32*[-0.88 -0.42 0.24];

    lmnVmsh = gseVinneredge*[L' M' N'];
    lmnVmsp = gseMSP*[L' M' N'];
     
    mspN = irf_norm(lmnVmsp);
    mshN = irf_norm(lmnVmsh);

    x = 90*[-1 1];
    y = 90*[-1 1];
    z = 90*[-1 1];
    
    funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
    funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
    funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);


    if exist('hmshN1'); delete(hmshN1); end
    if exist('hmshN2'); delete(hmshN2); end
    if exist('hmspN1'); delete(hmspN1); end
    if exist('hmspN2'); delete(hmspN2); end
    if exist('ht1'); delete(ht1); end
    if exist('ht2'); delete(ht2); end
    if exist('ht3'); delete(ht3); end
    if exist('ht4'); delete(ht4); end
   
    hold(hca,'on')
    if 1
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+70,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)-2.8,x*0,x,'k-');
    else
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+72,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)+20,x*0,x,'k-');
      hca.XLim = [0 50];
      hca.ZLim = [-10 40];
    end
    hold(hca,'off')

    ht1 = text(21.5,00,43,'MSH'); ht1.HorizontalAlignment = 'center'; ht1.FontSize = 13; ht1.Rotation = 55; 
    ht2 = text(2,0,50,'MSP'); ht2.HorizontalAlignment = 'center'; ht2.FontSize = 13; ht2.Rotation = -80;
    
    ht3 = text(16,0,-8,tintUTCstart(12:22)); ht3.HorizontalAlignment = 'center'; ht3.FontSize = 13;
    ht4 = text(16,0,72,tintUTCstop(12:22)); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
    
    hold(hca,'on') 
    quiver3(hca,17,0,4,0,0,5,2,'k')
    hold(hca,'off')
    ht4 = text(17,0,0,'time'); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
     
  end
end

%% Figure 2: Different comparison plot of more field components
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(13);
%irf_panel('delete1');
%irf_panel('delete2');
%irf_panel('delete3');
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 0 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end

if 1 % JeL
  hca = irf_panel('JeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).x,mvaJe2.tlim(tint).x,mvaJe3.tlim(tint).x,mvaJe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_e_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeM
  hca = irf_panel('JeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).y,mvaJe2.tlim(tint).y,mvaJe3.tlim(tint).y,mvaJe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_e_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeN
  hca = irf_panel('JeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).z,mvaJe2.tlim(tint).z,mvaJe3.tlim(tint).z,mvaJe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_e_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);  
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = EdJe?;',ic)
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'R_e\cdot J_e','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP  
  
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(mvaJe)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:9
  irf_legend(h(ii+3),legends{ii},[0.01 0.9],'color',[0 0 0])
end

%delete(h(1:3))
%%
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:35.00Z');
c_eval('h=irf_plot({mvaVExB?,mvaVe?perp,mvaVi?perp},''comp'');',ic)
irf_zoom(h(1:end),'x',tintZoom)
isub=1;
irf_legend(h(isub),{'v_{ExB}','v_{e,\perp}','v_{i,\perp}'},[0.02 0.95])
title(h(1),'Perpendicular velocities')
ylabel(h(1),'v_{L} (km/s)')
ylabel(h(2),'v_{M} (km/s)')
ylabel(h(3),'v_{N} (km/s)')


%%
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:35.00Z');
h = irf_plot(6);
if 1 % JeL
  hca = irf_panel('JeL');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaJ1.tlim(tint).x,mvaJ2.tlim(tint).x,mvaJ3.tlim(tint).x,mvaJ4.tlim(tint).x,mvaJcurl.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4','curl'},[0.02 0.95])
end
if 1 % JeM
  hca = irf_panel('JeM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaJ1.tlim(tint).y,mvaJ2.tlim(tint).y,mvaJ3.tlim(tint).y,mvaJ4.tlim(tint).y,mvaJcurl.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4','curl'},[0.02 0.95])
end
if 1 % JeN
  hca = irf_panel('JeN');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaJ1.tlim(tint).z,mvaJ2.tlim(tint).z,mvaJ3.tlim(tint).z,mvaJ4.tlim(tint).z,mvaJcurl.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4','curl'},[0.02 0.95])
end
if 1 % AvJL
  hca = irf_panel('AvJL');
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_plot(hca,{mvaAvJ.tlim(tint).x,mvaJcurl.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_legend(hca,{'fpi average','curl'},[0.02 0.95])
end
if 1 % AvJM
  hca = irf_panel('AvJM');
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_plot(hca,{mvaAvJ.tlim(tint).y,mvaJcurl.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_legend(hca,{'fpi average','curl'},[0.02 0.95])
end
if 1 % AvJN
  hca = irf_panel('AvJN');
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_plot(hca,{mvaAvJ.tlim(tint).z,mvaJcurl.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1b'))
  irf_legend(hca,{'fpi average','curl'},[0.02 0.95])
end

irf_zoom(h,'x',tintZoom)
title(h(1),'Current densities, from particle moments and magnetic field')

%% Supplemmentary figures: Figure MMS1, Particle distributions only perpendicular plane
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
ic = 4;
switch ic
  case 1
    tStart = EpochTT('2015-10-16T10:33:30.116567000Z');
    times = tStart + 0.03*[0:(nRows*nCols-1)];
    times = tStart + 0.03*[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
  case 2
    tStart = EpochTT('2015-10-16T10:33:29.838737000Z');
    times = tStart + 0.03*[0:(nRows*nCols-1)];
    times = tStart + 0.03*[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
  case 3
    tStart = EpochTT('2015-10-16T10:33:29.544423000Z');
    times = tStart + 0.03*[0:(nRows*nCols-1)];
    times = tStart + 0.03*[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
  case 4
    tStart = EpochTT('2015-10-16T10:33:29.548037000Z');
    times = tStart + 0.03*[0:(nRows*nCols-1)];
    times = tStart + 0.03*[0:2:31]+-1.5; 
end

nTimes = times.length;
% Initialize plots
h1 = irf_plot(10);
nRows = 4;
nCols = 4;
shiftRow = 2;
isub = 1;
clear h2;
for ii = (shiftRow*nCols+1):(nRows*(nCols+shiftRow)); h2(isub) = subplot(nRows+shiftRow,nCols,ii); isub = isub + 1;end

if 1
  hca = irf_panel(irf_ssub('B mms ?',ic));
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaB?.tlim(tintZoom).x,mvaB?.tlim(tintZoom).y,mvaB?.tlim(tintZoom).z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.02 0.9],'fontsize',12);
end
if 1
  hca = irf_panel(irf_ssub('Ve mms ?',ic));
  set(hca,'ColorOrder',mms_colors('axyz'))
  c_eval('irf_plot(hca,{mvaVe?par.tlim(tintZoom),mvaVe?perp.tlim(tintZoom).x,mvaVe?perp.tlim(tintZoom).y,mvaVe?perp.tlim(tintZoom).z},''comp'');',ic)
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('axyz'))
  irf_legend(hca,{'||','\perp_L','\perp_M','\perp_N'},[0.02 0.9],'fontsize',12);
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).e64.pitchangles(dmpaB?,24).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];  
  colormap(hca,'jet')
end

c_eval('dist = ePDist?;',ic)
c_eval('scpot = scPot?;',ic)
vlim = 12*1e3;
elevlim = 20;
strCMap = 'jet';
%energies =  [30 220];
projclim = [0 5];
palim = [1e-3 1e6];
skymapEnergy = [65 278];
hmark = irf_pl_mark(h1(1:2),times.epochUnix,'green');

isub = 1;
for ii = 1:nTimes
  time = times(ii);      

  % Get mean magnetic field direction
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);

  vectors = {hatB0,'B';hatE0,'E';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

  % Projection coordinate system
  if 1
    c_eval('z = gseB?.resample(time);',ic); 
    z = -z.data/norm(z.data); % B
    y = cross(z,N); % B x N
    x = cross(y,z); % N perp
    vlabels = {'B\times(N\times B)','B\times N','-B'};
    
    %z = -hatB0;
    %x = irf_norm(cross(hatB0,cross(N,hatB0)));
    %y = cross(z,x);   
    %vlabels = {'B\times(N\times B)','N_{\perp}','B'};
  elseif 1
    x = hatB0;
    y = hatExB0;
    z = cross(x,y);
  else
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
  end
  
  hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  mms.plot_projection(hca,dist.convertto('s^3/km^6'),'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'vlabel',vlabels,'scpot',scpot);
  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  hca.Title.String = '';
  colormap(hca,strCMap)
  hca.YLabel.String = '';
  hca.XLabel.String = '';
end


h2(1).YLabel.String = vlabels{2};
h2(end).XLabel.String = vlabels{1};

irf_zoom(h1(1:3),'x',tintZoom)
irf_zoom(h1(1),'y')
irf_plot_axis_align

for ii = 1:nRows*nCols
  originalPosition{ii} = h2(ii).Position;
end

hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
%hCB.Position = [h2(nCols).Position(1)+xWidth+0.02 h2(4).Position(2) 0.02 h2(4).Position(4)];

cmap = irf_colormap('space');
cmap = 'jet';
for ii = 1:nRows*nCols
  h2(ii).Position(3) = 0.15;originalPosition{ii}(3)*4;
  h2(ii).Position(4) = 0.15;originalPosition{ii}(4)*4;
  %h2(ii).Position(2) = originalPosition{ii}(2)+0.1;
  %h2(ii).Position(1) = originalPosition{ii}(1)+0.02;
  
  h2(ii).FontSize = 12;
  colormap(h2(11),cmap)
  h2(ii).XTickLabels = {}; 
  h2(ii).YTickLabels = {}; 
  h2(ii).XDir = 'reverse';
end

dx = 0.044;
for ii = [1:nCols:nTimes] % 1st column (left)
  %h2(ii).Position(2) = originalPosition{ii}(2)-0.08;
  h2(ii).Position(1) = originalPosition{ii}(1)+dx;
  h2(ii).YTick = [-10 0 10];
  h2(ii).YTickLabels = {'-10' '0' '10'};
end
for ii = [2:nCols:nTimes] % 2nd column
  %h2(ii).Position(2) = originalPosition{ii}(2)-0.08;
  h2(ii).Position(1) = originalPosition{ii}(1)+0;
end
for ii = [3:nCols:nTimes] % 3rd column
  %h2(ii).Position(2) = originalPosition{ii}(2)-0.08;
  h2(ii).Position(1) = originalPosition{ii}(1)-dx;
end
for ii = [4:nCols:nTimes] % 4th column
  %h2(ii).Position(2) = originalPosition{ii}(2)-0.08;
  h2(ii).Position(1) = originalPosition{ii}(1)-2*dx;
end

dy = 0.02;
for ii = 1:nCols % 1st row (top)
  h2(ii).Position(2) = originalPosition{ii}(2)-3*dy;
end
for ii = (nRows*1+1):(nRows*1+nCols) %
  h2(ii).Position(2) = originalPosition{ii}(2)-3*dy;
end
for ii = (nRows*2+1):(nRows*2+nCols) %
  h2(ii).Position(2) = originalPosition{ii}(2)-3*dy;
end
for ii = (nRows*3+1):(nRows*3+nCols) %
  h2(ii).Position(2) = originalPosition{ii}(2)-3*dy;
  h2(ii).XTick = [-10 0 10];
  h2(ii).XTickLabels = {'-10' '0' '10'};
end

%% Showing time matching
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
tintZoom = irf.tint('2015-10-16T10:33:28.00Z/2015-10-16T10:33:31.00Z');
dt = [0.00     -0.06     -0.18     -0.00];

h = irf_figure(6);

if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))    
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');  
  hca.YLabel.String = {'B_{M}','(nT)'};
end

if 1 % BM
  hca = irf_panel('BM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt);  
  hca.YLabel.String = {'B_{M}','(nT)'};
  irf_legend(hca,['dt = [' num2str(dt) ']'],[0.02 0.98])
end

if 1 % Ve perp L
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,L}','(km/s)'},'interpreter','tex');
end


if 1 % Ve perp L
  hca = irf_panel('Ve L  dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp','dt',dt);    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,L}','(km/s)'},'interpreter','tex');
  irf_legend(hca,['dt = [' num2str(dt) ']'],[0.02 0.98])
end

if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = 6*[-1 1];
end
if 1 % EN
  hca = irf_panel('EN dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp','dt',dt);
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,['dt = [' num2str(dt) ']'],[0.02 0.98])
  hca.YLim = 6*[-1 1];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(1:4),'y')
irf_plot_axis_align
h(1).Title.String = 'Time matching of magnetosheath side';
