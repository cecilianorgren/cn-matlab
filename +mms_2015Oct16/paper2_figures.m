%% Prepare data
units = irf_units;
%mms_2015Oct16.construct_particle_data % Make omni fluxes
ic = 1:4;
avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
avNi = (ni1+ni2.resample(ni1.time)+ni3.resample(ni1.time)+ni4.resample(ni1.time))/4; avNi.name = '<ni>';
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; gseAvB.name = '<B> (GSE)';
gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; gseAvVe.name = '<ve> (GSE)';
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.name = '<J> (GSE)';
vDe = gseGradPe.cross(gseAvB.resample(gseGradPe))/units.e/avNe/gseAvB.abs2*1e-9*1e-3; vDe.units = 'km/s'; vDe.name = 'v_De';
vDi = -gseGradPe.cross(gseAvB.resample(gseGradPe))/units.e/avNi/gseAvB.abs2*1e-9*1e-3; vDi.units = 'km/s'; vDi.name = 'v_Di';

% Calculate curvature of B
c_eval('gseMatR? = [gseR1brsttime.time.epochUnix gseR1brsttime.data];')
c_eval('gseMatB? = [gseB1.time.epochUnix gseB1.data];')
[gseCurvB,b]=c_4_grad('gseR?brsttime','gseB?','curvature');

% Plasma beta and magnetic pressure
c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

% Adiabatic invariant
c_eval('mu? = units.me*vte?perp*1e3*vte?perp*1e3/2/(gseB?.abs*1e-9)*1e9; mu?.units = ''nAm^2'';',ic)

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
c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)

c_eval('try [gseE?fastpar,gseE?fastperp] = irf_dec_parperp(gseB?.resample(gseE?fast),gseE?fast); end',ic); 
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?.resample(gseE?),gseE?);',ic)

c_eval('try gseVExB?fast = gseE?fast.cross(gseB?.resample(gseE?fast))/gseB?.resample(gseE?fast).abs2*1e3; gseVExB?fast.units = ''km/s''; end',ic);
c_eval('gseVExB? = gseE?.resample(gseB?).cross(gseB?)/gseB?.abs2*1e3; gseVExB?.units = ''km/s'';',ic)

%% Prepare pitchangle 
%c_eval('tic; ePitch? = ePDist?.pitchangles(dmpaB?,13); toc',ic)

%% Wave polarization analysis
% frequency = polarization.f;
% time = polarization.t;
% Bsum = polarization.bb_xxyyzzss(:,:,4);
% Esum = polarization.ee_xxyyzzss(:,:,4);
% Esum2D = polarization.ee_ss;
% ellipticity = polarization.ellipticity;
% dop = polarization.dop;
% thetak = polarization.k_tp(:,:,1);
% planarity = polarization.planarity;
% pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);

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
  case 5 % N: magnetosphere side normal derived from mms1 and mms4
    gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec34.time,M);
    gseNorm34 = gseVec34.cross(gseM);
    gseNormalMSP = gseNorm34/gseNorm34.abs;

    N = -gseNormalMSP.data;
    M = M;
    L = cross(M,N);
  case 6 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB4.tlim(irf.tint('2015-10-16T10:33:42.840Z/2015-10-16T10:33:53.517Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 7 % N: minimum variance of gradPe
    [out,l,v] = irf_minvar(gseGradPe.tlim(irf.tint('2015-10-16T10:33:47.331Z/2015-10-16T10:33:49.539Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 8 % N: minimum variance of Jcurl
    [out,l,v] = irf_minvar(gseJcurl.tlim(irf.tint('2015-10-16T10:33:44.317Z/2015-10-16T10:33:50.213Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 9 % N: minimum variance of curvB
    [out,l,v] = irf_minvar(gseCurvB.tlim(irf.tint('2015-10-16T10:33:44.240Z/2015-10-16T10:33:51.176Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 10 % N: minimum variance of '2nd (small) island' B1
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:47.181Z/2015-10-16T10:33:49.990Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
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
mvaCurvB = irf.ts_vec_xyz(gseCurvB.time,[gseCurvB.dot(L).data gseCurvB.dot(M).data gseCurvB.dot(N).data]); gseCurvB.units = gseCurvB.units;
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
c_eval('mvaJ?par = gseJ?par;')
c_eval('mvaJ?perp = irf.ts_vec_xyz(gseJ?perp.time,[gseJ?perp.dot(L).data gseJ?perp.dot(M).data gseJ?perp.dot(N).data]);')
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


%% Figure 1: 1 sc, B, vi, ve, def omni, pa, E & B wavelets
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:44.00Z/2015-10-16T10:33:51.00Z');
npanels = 9;
cmap = 'jet';
ic = 1;
boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);     
%tintZoom = tintZoom + [-5 5];
pdelete = [];
pshift = numel(pdelete);                      
h = irf_plot(npanels + pshift);
for ii = pdelete
    c_eval('irf_panel(''delete?'');',ii)
end 
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.e64.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,18).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_e_x','v_e_y','v_e_z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1200 700];  
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
if 0 % J moments 
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
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.tlim(tintZoom).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
  
  %hhleg=irf_legend(hca,irf_ssub('MMS ?',ic),[0.05 0.9],'color',[0 0 0],'fontsize',11);
  hcb.YLabel.String = {'Differential Energy Flux',hcb.YLabel.String{2}};
  hcb.YLabel.Position(2)=hcb.YLabel.Position(2)+0.5;
  hcb.YLabel.FontSize = 9;
  hca.CLim = [6 9];
  hca.YLim = [hca.YLim(1) 1e3];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 200];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 1000];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % E wavelet  
  hca=irf_panel('E wavelet'); %irf_plot(hca,dslE1brst);
  c_eval('[~,tind]=irf_tlim(wavE?.t,tintZoom.epochUnix'');',ic)
  %EpochTT(irf_time(wavE4.t,'epoch>utc'))
  c_eval('plotWavE? = wavE?; plotWavE?.t = plotWavE?.t(tind,1); plotWavE?.p{1} = plotWavE?.p{1}(tind,:);',ic)
  c_eval('irf_spectrogram(hca,plotWavE?,''log'',''donotfitcolorbarlabel'');',ic)
  %hcb(2) = colorbar('peer',hca);
  %hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic); irf_legend(hca,'f_{lh}',[0.95 0.2],'color','black');
  c_eval('hfce=irf_plot(hca,fce?,''white'');',ic); irf_legend(hca,'f_{ce}',[0.96 0.99],'color','white');
  % c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic); irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-8 -2]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
  colormap('jet');
end
if 1 % B wavelet
  hca=irf_panel('B wavelet'); 
  c_eval('[~,tind]=irf_tlim(wavB?.t,tintZoom.epochUnix'');',ic)
  c_eval('plotWavB? = wavB?; plotWavB?.t = plotWavB?.t(tind,1); plotWavB?.p{1} = plotWavB?.p{1}(tind,:);',ic)  
  c_eval('irf_spectrogram(hca,plotWavB?,''log'',''donotfitcolorbarlabel'');',ic)
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'E_y [(mV/m)^2/Hz]';
  hca.YScale = 'log';
  hold(hca,'on') 
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic); irf_legend(hca,'f_{lh}',[0.95 0.2],'color','black');
  c_eval('hfce=irf_plot(hca,fce?,''white'');',ic); irf_legend(hca,'f_{ce}',[0.96 0.99],'color','white');
  % c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic); irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');      
  hold(hca,'off')
  %legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end

if 1 % Ellipticity 
  hca=irf_panel('ellipticity'); 
  c_eval('[~,tind]=irf_tlim(polarization?.t,tintZoom.epochUnix'');',ic);  
  %c_eval('plotWavB? = wavB?; plotWavB?.t = plotWavB?.t(tind,1); plotWavB?.p{1} = plotWavB?.p{1}(tind,:);',ic)    
  c_eval('specrec=struct(''t'',polarization?.t(tind));',ic)
  c_eval('specrec.f=polarization?.f;',ic)
  c_eval('specrec.p=polarization?.ellipticity(tind,:);',ic)
  specrec.f_label='';
  specrec.p_label={'Ellipticity'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  hca.YScale = 'log';
  hold(hca,'on') 
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic); irf_legend(hca,'f_{lh}',[0.95 0.2],'color','black');
  c_eval('hfce=irf_plot(hca,fce?,''white'');',ic); irf_legend(hca,'f_{ce}',[0.96 0.99],'color','white');
  %c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic); irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  
  hold(hca,'off')  
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
  caxis(hca,[-1, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  %colormap(hca,irf_colormap('poynting'))
  colormap(hca,'jet')
end  
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end


irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:5]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

if 0% reset ylims 
hca = irf_panel('BN'); hca.YLim = [-7 7];
hca = irf_panel('BM'); hca.YTick = [-10 0 10];
hca = irf_panel('JM'); hca.YTick = [-1000 0 1000];
hca = irf_panel('E + vexB'); hca.YLim = [-4.5 6.5];
hca = irf_panel('EN'); hca.YLim = [-7.5 6.5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
%pshift = 3; % three deleted panels
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  %if labelsOutside
  %  h(ii+shift).Position(3) = h(ii+shift).Position(3)*0.88;
  %end  
  h(ii+pshift).YLabel.FontSize = 11;
end

for ii = []; h(ii).Visible = 'off'; end

h(3).YLim = [0 1e-9];
if 0 % Add labels for the different regions
  if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
  if exist('hleg_outflow','var'); delete(hleg_outflow); end
  if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
  hca = h(pshift+1);
  set(hca,'ColorOrder',mms_colors('11'))
  hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
  hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
  hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';
end
if 0 % Add borders/separtrices for the different regions
  hmark1 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
  hmark2 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
  for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; hmark1(ii).Color = [0.4 0.4 0.4]; end
  for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; hmark2(ii).Color = [0.4 0.4 0.4]; end
end
%% Figure 2: Plot figure with fields etc of 4 sc
figure(2)
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:44.00Z/2015-10-16T10:33:51.00Z');
tintZoom = irf.tint('2015-10-16T10:33:38.00Z/2015-10-16T10:33:53.00Z');
tint = tintZoom;
npanels = 12;
cmap = 'jet';
ic = 4;
boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);     
%tintZoom = tintZoom + [-5 5];
pshift = 0;                      
h = irf_plot(npanels + pshift);
for ii = 1:pshift
    c_eval('irf_panel(''delete?'');',ii)
end 
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
if 1 % Babs
  hca = irf_panel('Babs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
end
if 1 % curv B
  hca = irf_panel('curvB');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaCurvB.x,mvaCurvB.y,mvaCurvB.z},'comp');
  hca.YLabel.String = {'b \cdot \nabla b','(km^{-1})'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % 1/curv B
  hca = irf_panel('inv |curvB|');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{1/(mvaCurvB.abs)},'comp');
  hca.YLabel.String = {'|(b \cdot \nabla b)|^{-1}','(km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % JL
  hca = irf_panel('Jpar');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1par,mvaJ2par,mvaJ3par,mvaJ4par},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_{||}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % JL
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaJcurl.x},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % JN
  hca = irf_panel('JN');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.z,mvaJ2.z,mvaJ3.z,mvaJ4.z,mvaJcurl.z},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.tlim(tintZoom).specrec,''log'');',ic)  
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
  hca.CLim = [6 9];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all e');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).pitchangles(dmpaB?,20).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
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
if 0 % Ve par
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
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{mvaEVexB1.z,mvaEVexB2.z,mvaEVexB3.z,mvaEVexB4.z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  hca.YLim = [-10 10];
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


irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align
titleStr = ['L=[' num2str(L,'%5.2f') '], M=[' num2str(M,'%5.2f') '] N=[' num2str(N,'%5.2f') ']'];
h(1).Title.String = [titleStr ' (GSE)'];
if 0% reset ylims 
hca = irf_panel('BN'); hca.YLim = [-7 7];
hca = irf_panel('BM'); hca.YTick = [-10 0 10];
hca = irf_panel('JM'); hca.YTick = [-1000 0 1000];
hca = irf_panel('E + vexB'); hca.YLim = [-4.5 6.5];
hca = irf_panel('EN'); hca.YLim = [-7.5 6.5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

for ii = []; h(ii).Visible = 'off'; end

if 0 % Add labels for the different regions
  if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
  if exist('hleg_outflow','var'); delete(hleg_outflow); end
  if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
  hca = h(pshift+1);
  set(hca,'ColorOrder',mms_colors('11'))
  hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
  hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
  hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';
end
if 0 % Add borders/separtrices for the different regions
  hmark1 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
  hmark2 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
  for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; hmark1(ii).Color = [0.4 0.4 0.4]; end
  for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; hmark2(ii).Color = [0.4 0.4 0.4]; end
end
%% Figure 3: Adiabatic invariant
figure(3)
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:44.00Z/2015-10-16T10:33:51.00Z');
tintZoom = irf.tint('2015-10-16T10:33:35.00Z/2015-10-16T10:33:54.00Z');
tint = tintZoom;
npanels = 10;
cmap = 'jet';
ic = 4;
boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);     
%tintZoom = tintZoom + [-5 5];
pshift = 0;                      
h = irf_plot(npanels + pshift);
for ii = 1:pshift
    c_eval('irf_panel(''delete?'');',ii)
end

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.98 0.9],'fontsize',12);
end
if 1 % J par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  c_eval('irf_plot(hca,{mvaJ?perp.x.tlim(tint),mvaJ?perp.y.tlim(tint),mvaJ?perp.z.tlim(tint),mvaJ?par.tlim(tint)},''comp'');',ic)
  %hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'J_{L,\perp}','J_{M,\perp}','J_{N,\perp}','J_{||}'},[0.98 0.9],'fontsize',12);
  ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex');
end
if 1 % vte perp
  hca = irf_panel('veth perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vte?perp.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{te,\perp}','(km/s)'};    
end
if 1 % adiabatic invariant
  hca = irf_panel('mu');
  if 1 % 1 sc
    set(hca,'ColorOrder',mms_colors('1234'))
    c_eval('irf_plot(hca,{mu?.tlim(tint)},''comp'');',ic)
  else    
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{mu1.tlim(tint),mu2.tlim(tint),mu3.tlim(tint),mu4.tlim(tint)},'comp');    
  end
  hca.YLabel.String = {'\mu',['(' mu1.units ')']};
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
if 0 % Tepar/Teperp
  hca = irf_panel('Tepar/Teperp');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint)/((facTe?.yy+facTe?.zz)/2)},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}/T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 1 % beta
  hca = irf_panel('beta');
  if 1 % 1 sc
    set(hca,'ColorOrder',mms_colors('1234'))
    c_eval('irf_plot(hca,{beta?.tlim(tint)},''comp'');',ic)
  else % 4 sc    
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{beta1.tlim(tint),beta2.tlim(tint),beta3.tlim(tint),beta4.tlim(tint)},'comp');
  end
  hca.YLabel.String = {irf_ssub('\beta',ic)};
  hca.YScale = 'log';
  hca.YTick = [1e-1 1e0 1e1 1e2 1e3];
  hca.YLim = [1e-1 2e2];
  hca.MinorGridAlpha = 0;
end
if 1 % ne
  hca = irf_panel('ne');
  if 1 % 1 sc
    set(hca,'ColorOrder',mms_colors('1234'))
    c_eval('irf_plot(hca,{ne?.tlim(tint)},''comp'');',ic)
  else % 4 sc    
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  end
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % JL
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaJcurl.x},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 0 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 0 % JN
  hca = irf_panel('JN');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.z,mvaJ2.z,mvaJ3.z,mvaJ4.z,mvaJcurl.z},'comp');   
  % lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.tlim(tintZoom).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
  
  %hhleg=irf_legend(hca,irf_ssub('MMS ?',ic),[0.05 0.9],'color',[0 0 0],'fontsize',11);
  if 0
    hcb.YLabel.String = {'Differential Energy Flux',hcb.YLabel.String{2}};
    hcb.YLabel.Position(2)=hcb.YLabel.Position(2)+0.5;
    hcb.YLabel.FontSize = 9;
  end
  hca.CLim = [6 9];
  hca.YLim = [hca.YLim(1) 1e3];
end
if 1 % ePDist pa 32, low e 
  hca = irf_panel('e PA e32 deflux all low E 1');
  elim = [20 120];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 64, high e
  hca = irf_panel('e PA e32 deflux high E 2');  
  elim = [120 400];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa 32, low e 
  hca = irf_panel('e PA e32 deflux all low E 3');
  elim = [80 120];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa 64, high e
  hca = irf_panel('e PA e32 deflux high E 4');  
  elim = [120 400];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % Ve par
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


irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:5 7]),'y')
irf_plot_axis_align
titleStr = ['L=[' num2str(L,'%5.2f') '], M=[' num2str(M,'%5.2f') '] N=[' num2str(N,'%5.2f') ']'];
h(1).Title.String = [titleStr ' (GSE)'];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0;

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

for ii = []; h(ii).Visible = 'off'; end

if 0 % Add labels for the different regions
  if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
  if exist('hleg_outflow','var'); delete(hleg_outflow); end
  if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
  hca = h(pshift+1);
  set(hca,'ColorOrder',mms_colors('11'))
  hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
  hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
  hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';
end
if 0 % Add borders/separtrices for the different regions
  hmark1 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
  hmark2 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
  for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; hmark1(ii).Color = [0.4 0.4 0.4]; end
  for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; hmark2(ii).Color = [0.4 0.4 0.4]; end
end
%% Figure 4: Integrated spacecraft trajectory
% Quiver plot
figure(4)
h2 = subplot(1,1,1);
isub = 1;      
tintQuivers = irf.tint('2015-10-16T10:33:45.50Z/2015-10-16T10:33:50.00Z');
%tintQuivers = irf.tint('2015-10-16T10:33:44.00Z/2015-10-16T10:33:50.00Z');
%tintQuivers = irf.tint('2015-10-16T10:33:47.50Z/2015-10-16T10:33:50.00Z');
  
  times = mvaVe1.tlim(tintQuivers).time;  
  vel_selection = 2;
  clear timesVUTC gseVdata
  switch vel_selection
    case 1 % do velocities manually    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:42.50Z'; gseVdata(iv,:) = 37*[-0.14   -0.66   -0.74];
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:44.35Z'; gseVdata(iv,:) = 37*[-0.14   -0.66   -0.74];
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:45.55Z'; gseVdata(iv,:) = 37*[-0.14   -0.66   -0.74];
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:49.50Z'; gseVdata(iv,:) = 37*[-0.14   -0.66   -0.74];      
    case 2 % define velocities at certain times, and then interpolate to other times      
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:42.50Z'; gseVdata(iv,:) = 37*[-0.96 -0.26 -0.04]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:44.35Z'; gseVdata(iv,:) = 79*[0.41 0.73 -0.55]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:45.55Z'; gseVdata(iv,:) = 77*[0.29 0.13 -0.95]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:45.95Z'; gseVdata(iv,:) = 95*[0.40 0.43 -0.81]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:46.67Z'; gseVdata(iv,:) = 58*[0.58 0.46 -0.67]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:47.20Z'; gseVdata(iv,:) = 261*[-0.55 0.82 0.17]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:47.40Z'; gseVdata(iv,:) = 93*[0.28 0.96 -0.01]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:47.52Z'; gseVdata(iv,:) = 96*[0.18 0.97 0.18]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:47.63Z'; gseVdata(iv,:) = 144*[-0.15 0.93 -0.32]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:47.70Z'; gseVdata(iv,:) = 98*[0.46 0.85 -0.24]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:48.02Z'; gseVdata(iv,:) = 143*[-0.34 0.51 -0.79]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:48.35Z'; gseVdata(iv,:) = 123*[-0.78 0.43 -0.46]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:48.40Z'; gseVdata(iv,:) = 89*[-0.97 0.02 -0.24]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:48.48Z'; gseVdata(iv,:) = 78*[-0.99 0.13 -0.08]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:48.72Z'; gseVdata(iv,:) = 153*[-0.79 0.60 0.12]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:49.27Z'; gseVdata(iv,:) = 122*[-0.68 0.24 -0.69]; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:49.88Z'; gseVdata(iv,:) = 152*[0.15 0.65 -0.74]; 
  end
  gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
  gseV = gseV.resample(times);
  dt = times(2)-times(1);
  lmnV = gseV.data*[L' M' N'];
  
  
  c_eval('posR? = repmat(mvaRR?,times.length,1)-dt*cumsum(lmnV,1);')
  posR0 = mvaR0.resample(times)-mvaR0.resample(times.start).data; posR0 = posR0-dt*cumsum(lmnV,1); posR0=posR0.data;
  c_eval('posV? = mvaVe?.resample(times).data*0.6;')
  c_eval('posJ? = mvaJ?.resample(times).data;')
  c_eval('posB? = mvaB?.resample(times).data;')
  posCurvB = mvaCurvB.resample(times).data;
  posCurlJ = mvaJcurl.resample(times).data;
  posCurvBabs = mvaCurvB.resample(times).abs.data;
  posCurvB(posCurvBabs>0.010,:)=NaN;
  c_eval('posE? = mvaE?.resample(times).data;')
  c_eval('posRe? = mvaEVexB?.resample(times).data;')
  posGradPe = mvaGradPe.resample(times).data;
  
  hca = h2(isub); isub = isub +1;
  sclist = 1:4;
  hold(hca,'on')
  %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  %c_eval('plot_quivers(hca,[posV?(:,3) -posV?(:,2) posV?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  %c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
  %c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
  c_eval('plot_quivers(hca,[posB?(:,1) posB?(:,2) posB?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))',sclist)
  %c_eval('plot_quivers(hca,[posJ?(:,3) -posJ?(:,2) posJ?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  color = mms_colors('b');
  %plot_quivers(hca,[posCurvB(:,3) -posCurvB(:,2) posCurvB(:,1)],[posR0(:,3) -posR0(:,2) posR0(:,1)],color)
  %plot_quivers(hca,[posGradPe(:,3) -posGradPe(:,2) posGradPe(:,1)],[posR0(:,3) -posR0(:,2) posR0(:,1)],color)
  plot_quivers(hca,[posCurlJ(:,1) posCurlJ(:,2) posCurlJ(:,3)],[posR0(:,1) posR0(:,2) posR0(:,3)],color)
  %plot_quivers(hca,[posB1(:,3) -posB1(:,2) posB1(:,1)],[posR0(:,3) -posR0(:,2) posR0(:,1)],color)
  hold(hca,'off')

  hca.XLabel.String = 'L (km)';
  hca.YLabel.String = 'M (km)';
  hca.ZLabel.String = 'N (km)';
  
  %hca.ZDir = 'normal';
  %hca.YDir = 'normal';
  %hca.XDir = 'reverse';
  %axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  view(hca,[1 1 1])
  %hca.YLim = [-20 20];
  %axis(hca,'equal')
  axis(hca,'equal')
 % axis square
 
  if 0 % additional formatting
    hca.Box = 'on';
    hca.ZLim = [-10 70];
    hca.XLim = [-10 25];
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
  end
  if 0 % plot magnetosheath and magnetosphere boundary planes    
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
    hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+63,'k-.');
    hmspN1 = plot3(hca,funZ(x,y,mspN)-4,x*0,x,'k-');
    hold(hca,'off')

    ht1 = text(20,00,40,'MSH'); ht1.HorizontalAlignment = 'center'; ht1.FontSize = 13; ht1.Rotation = 55; 
    ht2 = text(-2,0,45,'MSP'); ht2.HorizontalAlignment = 'center'; ht2.FontSize = 13; ht2.Rotation = -80;
    
    ht3 = text(16,0,-8,tintUTCstart(12:22)); ht3.HorizontalAlignment = 'center'; ht3.FontSize = 13;
    ht4 = text(16,0,67,tintUTCstop(12:22)); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
    
    hold(hca,'on') 
    quiver3(hca,17,0,4,0,0,5,2,'k')
    hold(hca,'off')
    ht4 = text(17,0,0,'time'); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
     
  end
%% Figure 5: Plot sc positions and boundaires 
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
%% Figure 6: Detailed electron distributions where whislters are observed 
ic = 4;
nRows = 3;
nCols = 7;
nPanels = nRows*nCols;

isub = 0;
for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub + 1;         
    h(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
  end
end

tInd = {[1700 1701 1702 1703 1704]+6};
tInd = {[1690:5:1720],[1700:4:1720],[1700:4:1720],[1700:4:1760]+0*(177-60)};
c_eval('times = ePDist?.time;',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?;',ic)

isub=1;
for ii=1:nCols
  time = times(tInd{ic}(ii));
  timeUTC = time.utc;
  
  % Set up coordinate system for projection plots
  c_eval('x = dmpaB?.resample(times).resample(time); x = x.data/x.abs.data;',ic)
  y = cross(x,[1 0 0]); y = y/norm(y);
  z = cross(x,y);
  
  % Projection plot parameters
  vlim = 12*1e3;
  elevlim = 15;
  strCMap = 'jet';
  projclim = [-1 5];
  
  % Projections
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
  %hca.Title.String = '';
  hca.Title.String = timeUTC(12:23);
  colormap(hca,strCMap)                
  hca.XLabel.String = 'v_{\perp,1}';
  hca.YLabel.String = 'v_{\perp,2}';
  
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
  hca.Title.String = '';
  %hca.Title.String = timeUTC(12:23);
  colormap(hca,strCMap)                
  hca.XLabel.String = 'v_{\perp,2}';
  hca.YLabel.String = 'v_{||}';
  
  % Pitch angle distributions
  hca = h(isub); isub = isub + 1; 
    % Add a reference maxwellian, so one can see the evolution easier
  n = [20];
  m = [0];
  t = [35];
  vd = [0];
  d = [1];
  a1 = [1];
  a2 = [0];
  toPlot = [1];
  Ee = ePitch1.depend{1}(tInd{ic}(ii),:);
  ve = cn_eV2v(Ee,'eV');
  f = cn.maxwellian(ve,35,20,0,'e','3D'); % f = cn.maxwellian(v,T,n,vd,species,optional);  
  
  loglog(hca,Ee,f*1e18,'k:'); % 1e18: m^-6 ->km^-6
  %hh=whamp.plot_f(hca,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
  
  hold(hca,'on')
  unitscale = 1e30; % cm^-6 ->km^-6
  plot(hca,ePitch1.depend{1}(tInd{ic}(ii),:),ePitch1.data(tInd{ic}(ii),:,1)*unitscale,...
           ePitch1.depend{1}(tInd{ic}(ii),:),ePitch1.data(tInd{ic}(ii),:,7)*unitscale,...
           ePitch1.depend{1}(tInd{ic}(ii),:),ePitch1.data(tInd{ic}(ii),:,13)*unitscale)
  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLim = [1e-32 2e-25]*unitscale;
  
  hca.YLabel.String = 'f_e (s^3 km^{-6})';
  hca.XTick = [1e1 1e2 1e3 1e4];
  hleg = irf_legend(hca,{'0';'90';'180'},[0.98 0.98]);

  hca.Title.String = '';
end

hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.04;

h(1).Title.String = {irf_ssub('MMS 1',ic),h(1).Title.String};

