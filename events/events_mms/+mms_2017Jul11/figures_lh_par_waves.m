tint_zoom1 = irf.tint('2017-07-11T22:33:17.00Z/2017-07-11T22:33:25.00Z'); % wher
tint_zoom2 = irf.tint('2017-07-11T22:33:19.00Z/2017-07-11T22:33:21.00Z'); %

%% Rotate into LMN coordinates
ic = 1:4;
% Set up coordinate system
% Default
%[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:19.611Z/2015-11-12T07:19:22.641Z')));
%L = v(1,:); M = v(2,:); N = v(3,:);
%coordLabels = {'L','M','N'};
%lmn = [N;-M;L];

% [out,l,v] = irf_minvar(gseAvB.tlim(irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z')));
% L = v(1,:); M = v(2,:); N = v(3,:);
% coordLabels = {'L','M','N'};
% lmn = [L;M;N];

% dt_timing = [0.0000   0.0141   0.0104   0.0078];
% v_timing = 1.34e+03*[-0.85 -0.48  0.24];
% c_eval('tsV?_timing = irf.ts_vec_xyz(gseE?.time,repmat(v_timing,gseE?.length,1));',ic)
% v_direction = irf_norm(v_timing);
% v_amplitude = sqrt(sum(v_timing.^2));
% 
dt_timing = [0.0000    0.0247    0.0216    0.0093];
v_timing = 737*[-0.67 -0.53  0.51];
v_direction = irf_norm(v_timing);
% tintLH = irf.tint('2017-07-11T22:33:28.60Z/2017-07-11T22:33:29.30Z'); % LH wave packet
newz = irf_norm(mean(gseB1.tlim(tint_zoom2).data,1)); % B
newx = cross(newz,cross(v_direction,newz)); % J direction
newy = cross(newz,newx); % normal to v/j and B
newxyz = [newx;newy;newz];

lmn = newxyz;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);

tintFluxrope = irf.tint('2017-07-11T22:33:11.557111083Z/2017-07-11T22:33:41.157944824Z');
%[~,lB,lmnB]=irf_minvar(gseB1);
%lmn = lmnB;

if 0
N = newy; % direction of Eperp
M = cross(N,cross(lmnB(2,:),N));
L = cross(M,N);
lmn = [L;M;N];
end

disp(sprintf('L = [%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L,M,N))
% Rotate data
c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';')
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';')
c_eval('mvaB?scm = gseB?scm*lmn''; mvaB?scm.name = ''B LMN'';')
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
c_eval('mvaRR? = mvaR?-mvaR0; mvaRR? = mvaRR?.resample(irf_time(''2015-11-12T07:19:21.000Z'',''utc>epochTT'')).data;',ic)
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

%% Reduced dist
%eDist = ePDist1.tlim(tint);
ic = 1;
c_eval('eDist = ePDist?.tlim(tint_zoom1);',ic)

if 0 % remove background
  nSecondary = [5];
  nPhoto = 0;
  %[eDist_nobg] = mms.remove_edist_background(eDist_orig);
  c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''Nphotoe_art'',nPhoto,''ZeroNaN'',0);',1:numel(nSecondary))
end



vg = (-100:2:100)*1e3;
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist).norm;
zhat = irf.ts_vec_xyz(ePara.time,repmat(M,ePara.length,1));
ePerp1 = zhat.cross(ePara).norm;
ePerp2 = ePara.cross(ePerp1).norm;

%%
eint = [000 40000];
vint = [-Inf Inf];
scpot = scPot1.resample(eDist);
lowerelim = 0;
nMC = 500;

tic; ef1D_para = eDist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp1 = eDist.reduce('1D',ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_perp2 = eDist.reduce('1D',ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); toc % reduced distribution along B
%tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
%tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
%tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc


%% Figure
npanels = 9;
h = irf_plot(npanels);
nrows = 2;
ncols = 2;
%[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.5,'vertical');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z},'comp');
  hca.YLabel.String = {'B (GSE)','(nT)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % Epar 4sc
  hca = irf_panel('Epar 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % Eperp x 4sc
  hca = irf_panel('Eperp x 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.x,gseE2perp.x,gseE3perp.x,gseE4perp.x},'comp');
  hca.YLabel.String = {'E_{\perp,x}','(mV/m)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % Eperp y 4sc
  hca = irf_panel('Eperp y 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.y,gseE2perp.y,gseE3perp.y,gseE4perp.y},'comp');
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % Eperp z 4sc
  hca = irf_panel('Eperp z 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.z,gseE2perp.z,gseE3perp.z,gseE4perp.z},'comp');
  hca.YLabel.String = {'E_{\perp,z}','(mV/m)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % vepar 4sc
  hca = irf_panel('ve par 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par,gseVe2par,gseVe3par,gseVe4par},'comp');
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end

if 1 % e psd vpara
  hca = irf_panel('fred para');
  [~,hcb] = irf_spectrogram(hca,ef1D_para.specrec('velocity_1D','10^3 km/s'));
  if 0 % electron moment along projection line
    hold(hca,'on')
    irf_plot(hca,{lineVe*1e-3},'comp')
    %irf_plot(hca,gseVi1)  
    hold(hca,'off')
  end
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  hcb.YLabel.String = {'log_{10}f(v_{||})','(s/m^4)'};
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.1],'color',1*[1 1 1])
  if 0
    if unique(ef1D.ancillary.lowerelim) == 1
      irf_legend(hca,['E_{e} >' num2str(unique(ef1D.ancillary.lowerelim)) 'V_{sc}'],[0.99 0.1],'color',1*[1 1 1])
    end
  end
end
if 1 % e psd vperp1
  hca = irf_panel('fred perp1');
  [~,hcb] = irf_spectrogram(hca,ef1D_perp1.specrec('velocity_1D','10^3 km/s'));
  if 0 % electron moment along projection line
    hold(hca,'on')
    irf_plot(hca,{lineVe*1e-3},'comp')
    %irf_plot(hca,gseVi1)  
    hold(hca,'off')
  end
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,\perp,1}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  hcb.YLabel.String = {'log_{10}f(v_{\perp,1})','(s/m^4)'};
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.1],'color',1*[1 1 1])
  if 0
    if unique(ef1D.ancillary.lowerelim) == 1
      irf_legend(hca,['E_{e} >' num2str(unique(ef1D.ancillary.lowerelim)) 'V_{sc}'],[0.99 0.1],'color',1*[1 1 1])
    end
  end
end
if 1 % e psd vperp2
  hca = irf_panel('fred perp2');
  [~,hcb] = irf_spectrogram(hca,ef1D_perp2.specrec('velocity_1D','10^3 km/s'));
  if 0 % electron moment along projection line
    hold(hca,'on')
    irf_plot(hca,{lineVe*1e-3},'comp')
    %irf_plot(hca,gseVi1)  
    hold(hca,'off')
  end
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,perp,2}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  hcb.YLabel.String = {'log_{10}f(v_{perp,1})','(s/m^4)'};
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.1],'color',1*[1 1 1])
  if 0
    if unique(ef1D.ancillary.lowerelim) == 1
      irf_legend(hca,['E_{e} >' num2str(unique(ef1D.ancillary.lowerelim)) 'V_{sc}'],[0.99 0.1],'color',1*[1 1 1])
    end
  end
end
if 0 % ePDist deflux omni  
  hca = irf_panel('e DEF omni');  
  irf_spectrogram(hca,ePDist1.tlim(tint_zoom1).deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  %irf_plot(hca,gseTe1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end

irf_zoom(h,'x',tint_zoom1)
irf_zoom(h,'y')
irf_plot_axis_align

fred_vlim = 60*[-1 1];
hca = irf_panel('fred para'); hca.YLim = fred_vlim;
hca = irf_panel('fred perp1'); hca.YLim = fred_vlim;
hca = irf_panel('fred perp2'); hca.YLim = fred_vlim;

%% Figure, hodograms of mvaVe, 4 panels
rows = 2;
cols = 2;
panels = nrows*ncols;
h = setup_subplots(rows,ncols);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  v1 = mvaVe1.tlim(tint_zoom2).x.data;
  v2 = mvaVe1.tlim(tint_zoom2).y.data;
  plot(hca,v1,v2)  
end  
if 1
  hca = h(isub); isub = isub + 1;
  v1 = mvaVe2.tlim(tint_zoom2).x.data;
  v2 = mvaVe2.tlim(tint_zoom2).y.data;
  plot(hca,v1,v2)  
end  
if 1
  hca = h(isub); isub = isub + 1;
  v1 = mvaVe3.tlim(tint_zoom2).x.data;
  v2 = mvaVe3.tlim(tint_zoom2).y.data;
  plot(hca,v1,v2)  
end  
if 1
  hca = h(isub); isub = isub + 1;
  v1 = mvaVe4.tlim(tint_zoom2).x.data;
  v2 = mvaVe4.tlim(tint_zoom2).y.data;
  plot(hca,v1,v2)  
end  
%% Figure, hodograms of mvaVe, 1 panels
nrows = 1;
ncols = 3;
panels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

utc_tint1 = tint_zoom2(1).utc;
utc_tint2 = tint_zoom2(1).utc;
utc_tint = tint_zoom2.utc;
utc_tint = utc_tint(:,1:19);

vlim = 4999*[-1 1];
if 1 % xy
  hca = h(isub); isub = isub + 1;
  vx1 = mvaVe1.tlim(tint_zoom2).x.data;
  vy1 = mvaVe1.tlim(tint_zoom2).y.data;
  vx2 = mvaVe2.tlim(tint_zoom2).x.data;
  vy2 = mvaVe2.tlim(tint_zoom2).y.data;
  vx3 = mvaVe3.tlim(tint_zoom2).x.data;
  vy3 = mvaVe3.tlim(tint_zoom2).y.data;
  vx4 = mvaVe4.tlim(tint_zoom2).x.data;
  vy4 = mvaVe4.tlim(tint_zoom2).y.data;
  plot(hca,vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_{\perp,1} \sim wave propagationg direction';
  hca.YLabel.String = 'v_{\perp,2} \sim boundary normal direction';
end  
if 1 % yz
  hca = h(isub); isub = isub + 1;
  vx1 = mvaVe1.tlim(tint_zoom2).y.data;
  vy1 = mvaVe1.tlim(tint_zoom2).z.data;
  vx2 = mvaVe2.tlim(tint_zoom2).y.data;
  vy2 = mvaVe2.tlim(tint_zoom2).z.data;
  vx3 = mvaVe3.tlim(tint_zoom2).y.data;
  vy3 = mvaVe3.tlim(tint_zoom2).z.data;
  vx4 = mvaVe4.tlim(tint_zoom2).y.data;
  vy4 = mvaVe4.tlim(tint_zoom2).z.data;
  plot(hca,vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_{\perp,2} \sim boundary normal direction';
  hca.YLabel.String = 'v_{||}';
end  
if 1 % xz
  hca = h(isub); isub = isub + 1;
  vx1 = mvaVe1.tlim(tint_zoom2).x.data;
  vy1 = mvaVe1.tlim(tint_zoom2).z.data;
  vx2 = mvaVe2.tlim(tint_zoom2).x.data;
  vy2 = mvaVe2.tlim(tint_zoom2).z.data;
  vx3 = mvaVe3.tlim(tint_zoom2).x.data;
  vy3 = mvaVe3.tlim(tint_zoom2).z.data;
  vx4 = mvaVe4.tlim(tint_zoom2).x.data;
  vy4 = mvaVe4.tlim(tint_zoom2).z.data;
  plot(hca,vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_{\perp,1} \sim wave propagationg direction';
  hca.YLabel.String = 'v_{||}';
end

h(2).Title.String = [utc_tint(1,:) ' - ' utc_tint(2,:)];
%% Quivers along trajectory

%% mms.lhwaveanalysis
Tint = irf.tint('2017-07-11T22:33:28.70Z/2017-07-11T22:33:29.30Z'); % dont know
Exyz = gseE1;
Bscm = gseB1scm;
Bxyz = gseB1;
ne = ne1;
[phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(Tint,Exyz,Bscm,Bxyz,ne,'plot',1);

%% Figure
npanels = 7;
nrows = 2;
ncols = 2;
%h = irf_plot(npanels);
[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.5,'vertical');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z},'comp');
  hca.YLabel.String = {'B (GSE)','(nT)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % BDC and BAC
  hca = irf_panel('B dc ac abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.abs,dcB1.abs},'comp');
  hca.YLabel.String = {'|B| (GSE)','(nT)'};
  irf_legend(hca,{'',sprintf('f<%g',ffilt),'z'},[0.98 0.98])
end
if 0 % gseE
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.x,gseE1.y,gseE1.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % gseE, lowpass
  hca = irf_panel('E lowpass');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{plotE.x,plotE.y,plotE.z},'comp');
  hca.YLabel.String = {'E (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 1 % gseVixB
  hca = irf_panel('vixB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVixB1.x,-gseVixB1.y,-gseVixB1.z},'comp');
  hca.YLabel.String = {'-v_i\times B (GSE)','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % gseEperp+vixB
  hca = irf_panel('E+vixB');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseEperpVixB1.x,gseEperpVixB1.y,gseEperpVixB1.z},'comp');
  hca.YLabel.String = {'E_\perp+v_i\times B','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 0 % gseE, lowpass, vvixB x
  hca = irf_panel('E vixb x');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.x,plotE.x,-gseVixB1.x},'comp');
  hca.YLabel.String = {'E_x (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_x','E_x','-v_ixB_x'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 0 % gseE, lowpass, vvixB y
  hca = irf_panel('E vixb y');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.y,plotE.y,-gseVixB1.y},'comp');
  hca.YLabel.String = {'E_y (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_y','E_y','-v_ixB_y'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 0 % gseE, lowpass, vvixB z
  hca = irf_panel('E vixb z');
  set(hca,'ColorOrder',mms_colors('1234'))
  fhigh = 0.5;
  plotE = gseE1.filt(0,fhigh,[],3);
  irf_plot(hca,{-gseVexB1.z,plotE.z,-gseVixB1.z},'comp');
  hca.YLabel.String = {'E_z (GSE)','(mV/m)'};
  irf_legend(hca,{'-v_exB_z','E_z','-v_ixB_z'},[0.98 0.98])
  irf_legend(hca,sprintf('f<%g Hz',fhigh),[0.02 0.98],'color',[0 0 0])
end
if 1 % iPDist deflux omni  
  hca = irf_panel('i DEF omni');  
  irf_spectrogram(hca,iPDist1.deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 0 % Ti part
  %%
  figure;
  
  hca = irf_panel('i Ti part');  
  irf_spectrogram(hca,Ti1partdiff.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  Tsum = irf.ts_scalar(Ti1part.time,cumsum(Ti1part.data,2));
  irf_plot(hca,{gseTi1.trace/3,Tsum},'comp');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 0 % iPDist dpflux omni feeps 
  hca = irf_panel('i PEF omni feeps');  
  irf_spectrogram(hca,feepsi1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 0 % iPDist dpflux omni EIS
  hca = irf_panel('i PEF EIS omni');  
  irf_spectrogram(hca,eisi1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};  
end
if 0 % iPDist dpflux omni  FPI
  hca = irf_panel('i PEF omni');  
  irf_spectrogram(hca,iPDist1.dpflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTi1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_i^{FPI}','(eV)'};  
end
if 0 % iPDist dpflux omni, HPCA  
  hca = irf_panel('i DPF omni HPCA');
  irf_spectrogram(hca,hplusOmni1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,Thplus1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_{H+}^{HPCA}','(eV)'}; 
  irf_legend(hca,{'HCPA'},[0.02,0.02],'color',[0 0 0])
end
if 0 % oPDist deflux omni, HPCA  
  hca = irf_panel('o DPF omni HPCA');  
  irf_spectrogram(hca,oplusOmni1.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,Toplus1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_{O+}','(eV)'}; 
  irf_legend(hca,{'HCPA'},[0.02,0.02],'color',[0 0 0])
end
if 1 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  irf_spectrogram(hca,ePDist1.deflux.omni.specrec,'log');
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  irf_plot(hca,gseTe1.trace/3,'k');
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 1 % Ti, Te
  hca = irf_panel('T');
  set(hca,'ColorOrder',mms_colors('2134'))
  irf_plot(hca,{gseTi1.trace/3,Thplus1.trace/3,gseTe1.trace/3},'comp');
  hca.YLabel.String = {'T','(eV)'};
  irf_legend(hca,{'fpi i','hpca h+','e'},[0.98 0.98])
end
if 0 % ni,ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('2134'))
  irf_plot(hca,{ni1,nhplus1,noplus1,ne1},'comp');
  hca.YLabel.String = {'n','(cm^{-3})'};
  irf_legend(hca,{'fpi i','hpca h+','hpca o+','fpi e'},[0.98 0.98])
end

irf_plot_axis_align
irf_pl_mark(h,tint_ti)
irf_zoom(h,'x',tint)

% ion
omni_fpi.e = iPDist1.tlim(tint_ti).depend{1}(1,:);
omni_fpi.d = nanmean(iPDist1.tlim(tint_ti).dpflux.omni.data,1);
omni_hpca.e = hplusOmni1.depend{1}(:);
omni_hpca.d = nanmean(hplusOmni1.tlim(tint_ti).data,1);
%omni_eis.e = Omnifluxion_epd_eis_brst_l2.depend{1}(1,:);
%omni_eis.d = nanmean(Omnifluxion_epd_eis_brst_l2.data,1)*1e-3; % 1/keV -> /eV
omni_feeps.e = feepsi1.tlim(tint_ti).depend{1}(:);
omni_feeps.d = nanmean(feepsi1.tlim(tint_ti).data,1)*1e-3; % 1/keV -> /eV

ti_av = mean(gseTi1.tlim(tint_ti).trace.data/3,1);
vi_av = mean(gseVi1.tlim(tint_ti).abs.data,1);
Evi_av = units.mp*(vi_av*1e3).^2/2/units.eV;

% electrons
omnie_fpi.e = ePDist1.tlim(tint_ti).depend{1}(1,:);
omnie_fpi.d = nanmean(ePDist1.tlim(tint_ti).dpflux.omni.data,1);
%omni_eis.e = Omnifluxion_epd_eis_brst_l2.depend{1}(1,:);
%omni_eis.d = nanmean(Omnifluxion_epd_eis_brst_l2.data,1)*1e-3; % 1/keV -> /eV
%omni_feeps.e = feepsi1.tlim(tint_ti).depend{1}(:);
%omni_feeps.d = nanmean(feepsi1.tlim(tint_ti).data,1)*1e-3; % 1/keV -> /eV

te_av = mean(gseTe1.tlim(tint_ti).trace.data/3,1);
ve_av = mean(gseVe1.tlim(tint_ti).abs.data,1);
Eve_av = units.me*(ve_av*1e3).^2/2/units.eV;

