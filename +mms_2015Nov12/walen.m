% Wal?n test
units = irf_units;
%c_eval('facTe?fast = mms.rotate_tensor(gseTe?fast,''fac'',gseB?srvy); facTe?fast,name = ''fac Te fast'';',ic)

%% Magnetosheat reference point.
timeMSH = irf_time('2015-11-12T07:25:00.000Z','utc>epochTT');
timeMSH = irf_time('2015-11-12T06:32:00.000Z','utc>epochTT');
dt = 60*3;
tintMSH = irf.tint('2015-11-12T06:44:00.00Z',dt);
tintMSH = irf.tint('2015-11-12T07:24:40.00Z',dt);
tintMSH = irf.tint('2015-11-12T07:28:80.00Z',dt);
meantime = tintMSH.start+dt/2;
c_eval('dP? = irf.ts_scalar(facPe?fast.time,facPe?fast.xx.data-0.5*(facPe?fast.yy.data+facPe?fast.zz.data));',ic)
c_eval('Alfa? = dP?/units.mu0/gseB?srvy.resample(dP?).abs2;',ic)
c_eval('Ve? = gseVe?fast;',ic)
c_eval('Vi? = gseVi?fast;',ic)
c_eval('B? = gseB?srvy;',ic)

if 0
  timeused = timeMSH;
  c_eval('mshAlfa? = Alfa?.resample(timeMSH);',ic)
  c_eval('mshVe? = gseVe?fast.resample(timeMSH);',ic)
  c_eval('mshVi? = gseVi?fast.resample(timeMSH);',ic)
  c_eval('mshB? = gseB?srvy.resample(timeMSH);',ic)
else
  timeused = tintMSH;
  c_eval('mshAlfa? = Alfa?.tlim(tintMSH); mshAlfa? = irf.ts_scalar(meantime,mean(mshAlfa?.data,1));',ic)
  c_eval('mshVe? = gseVe?fast.tlim(tintMSH);  mshVe? = irf.ts_vec_xyz(meantime,mean(mshVe?.data,1));',ic)
  c_eval('mshVi? = gseVi?fast.tlim(tintMSH);  mshVi? = irf.ts_vec_xyz(meantime,mean(mshVi?.data,1));',ic)
  c_eval('mshB? = gseB?srvy.tlim(tintMSH); mshB? = irf.ts_vec_xyz(meantime,mean(mshB?.data,1));',ic)  
end

c_eval('dVe?obs = Ve?-mshVe?.data; dVe?obs.name = ''delta ve obs'';',ic)
c_eval('dVe?th = (gseB?srvy.resample(Alfa?)*(-Alfa?+1)-mshB?.data*(1-mshAlfa?.data))/ne?fast.sqrt/units.mu0*1e-6; dVe?th.name = ''delta ve th'';',ic)
c_eval('dVi?obs = Vi?-mshVi?.data; dVi?obs.name = ''delta vi obs'';',ic)
c_eval('dVi?th = (gseB?srvy.resample(Alfa?)*(-Alfa?+1)-mshB?.data*(1-mshAlfa?.data))/ni?fast.sqrt/units.mu0*1e-6; dVi?th.name = ''delta vi th'';',ic)

c_eval('walAngle?e = irf.ts_scalar(dVe?obs.time,acosd(dVe?obs.dot(dVe?th.resample(dVe?obs)/dVe?th.abs/dVe?obs.abs).data)); walAngle?.name = ''Walen angle e'';',ic)
c_eval('walR?e = dVe?obs.abs/dVe?th.resample(dVe?obs).abs;',ic) % walR.name = ''Walen ratio'';

c_eval('walAngle?i = irf.ts_scalar(dVi?obs.time,acosd(dVi?obs.dot(dVi?th.resample(dVi?obs)/dVi?th.abs/dVi?obs.abs).data)); walAngle?i.name = ''Walen angle i'';',ic)
c_eval('walR?i = dVi?obs.abs/dVi?th.resample(dVi?obs).abs;',ic) % walR.name = ''Walen ratio'';

%h = irf_plot({B1,Ve1,gseVi1fast,dV1obs,dV1th,walAngle1e,walR1e,walAngle1i,walR1i});
%mshRefMark = irf_pl_mark(h,timeused.epochUnix');

%% Make plot, fast data
tint = irf.tint('2015-11-12T06:05:00.00Z/2015-11-12T07:55:00.00Z');
%tint = irf.tint('2015-11-12T06:55:00.00Z/2015-11-12T07:32:00.00Z');
npanels = 7;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?srvy.x.tlim(tint),gseB?srvy.y.tlim(tint),gseB?srvy.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?fast},''comp'');',ic)
  hca.YLabel.String = {'n_i','(cm^{-3})'}; 
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?fast.x.tlim(tint),gseVi?fast.y.tlim(tint),gseVi?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?fast.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?fast.deflux.tlim(tint).pitchangles(dmpaB?srvy,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 1 % Ti par perp
  hca = irf_panel('Ti');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTi?fast.xx.tlim(tint),(facTi?fast.yy+facTi?fast.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];  
end
if 1 % Walen ratio
  hca = irf_panel('Walen angle');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walAngle?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\theta_W'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [0 180];  
end
if 1 % Walen angle
  hca = irf_panel('Walen ratio');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walR?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'R_W'};
  set(hca,'ColorOrder',mms_colors('1'))
  hca.YLim = [0 5];
end
if 0 % alfa
  hca = irf_panel('alpha');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Alfa?.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\alpha'};
  set(hca,'ColorOrder',mms_colors('1'))  
end


irf_zoom(h,'x',tint)
irf_zoom(h([1:3 5]),'y')
irf_plot_axis_align
mshRefMark = irf_pl_mark(h,timeused.epochUnix');
h(1).Title.String = irf_ssub('MMS ?',ic);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Make plot, fast data
tint = irf.tint('2015-11-12T06:05:00.00Z/2015-11-12T07:55:00.00Z');
%tint = irf.tint('2015-11-12T06:55:00.00Z/2015-11-12T07:32:00.00Z');
npanels = 8;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?srvy.x.tlim(tint),gseB?srvy.y.tlim(tint),gseB?srvy.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?fast,ne?fast},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'}; 
  irf_legend(hca,{'n_i','n_e'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?fast.x.tlim(tint),gseVi?fast.y.tlim(tint),gseVi?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?fast.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?fast.x.tlim(tint),gseVe?fast.y.tlim(tint),gseVe?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-300 220];
end
if 1 % ePDist omni 64
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?fast.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?fast.deflux.tlim(tint).pitchangles(dmpaB?srvy,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 1 % Ti par perp
  hca = irf_panel('Ti');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTi?fast.xx.tlim(tint),(facTi?fast.yy+facTi?fast.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];  
end
if 0 % Walen ratio
  hca = irf_panel('Walen angle');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walAngle?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\theta_W'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [0 180];  
end
if 0 % Walen angle
  hca = irf_panel('Walen ratio');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walR?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'R_W'};
  set(hca,'ColorOrder',mms_colors('1'))
  hca.YLim = [0 5];
end
if 0 % alfa
  hca = irf_panel('alpha');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Alfa?.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\alpha'};
  set(hca,'ColorOrder',mms_colors('1'))  
end


irf_zoom(h,'x',tint)
irf_zoom(h([1:3 5 7]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Wenyas method
%% Example code to perform Walen test; only for burst mode MMS data.
%   Reference: Phan et al. [AG, 2004]
%%  1. basic
    ic = 1;
    Jsign = -1;         % 1/-1 for jet direction;    
    %time = irf_time('2015-11-12T06:43:00.000Z','utc>epochtt');
    vec1 = [1, 0, 0];       vec2 = [0, 1, 0];       vec3 = [0, 0, 1];       % in GSE
    timeinterval = 2;
    switch timeinterval
      case 1
        tint = irf.tint('2015-11-12T06:29:00.00Z/2015-11-12T06:35:00.00Z');
        tint_ref = irf.tint('2015-11-12T06:30:00.00Z/2015-11-12T06:31:00.00Z');
        tint_walen = irf.tint('2015-11-12T06:32:00.00Z/2015-11-12T06:35:00.00Z');
      case 2
        tint = irf.tint('2015-11-12T06:42:05.00Z/2015-11-12T06:43:14.00Z');
        tint_ref = irf.tint('2015-11-12T06:43:00.00Z/2015-11-12T06:43:10.00Z');
        tint_walen = irf.tint('2015-11-12T06:42:45.00Z/2015-11-12T06:42:55.00Z');
        vec1 = [0.3178    0.0351    0.9475];      
        vec2 = [-0.2694   -0.9548    0.1258];       
        vec3 = [0.9091   -0.2952   -0.2940];       % in GSE               
      otherwise        
        tint = irf.tint('2015-11-12T06:42:05.000000Z/2015-11-12T06:43:10.000000Z');             % plot
        tint_ref = irf.tint('2015-11-12T07:21:30.000000Z/2015-11-12T07:24:00.000000Z');         % reference region     
        tint_walen = irf.tint('2015-11-12T07:19:00.000000Z/2015-11-12T07:21:30.000000Z');       % test region
    end
    unit_v2w = 1e6 * 0.5 * 1.673e-27 / 1.6e-19;
    time = tint(1)+10;
%%  2. load brst data, choose brst or fast
    % 2.1. PSD
    c_eval('ifile? = mms.get_filepath([''mms?_fpi_brst_l2_dis-dist''], time);', ic);
    c_eval('[iPDist?, iPDistError?] = mms.make_pdist(ifile?);', ic); 
    % 2.2. moments
    c_eval('Tint = irf.tint(iPDist?.time.start, iPDist?.time.stop);', ic);
    c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',Tint, ?);',ic);
    c_eval('Vi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint,?);',ic);
    c_eval('Pi? = mms.get_data(''Pi_dbcs_fpi_brst_l2'',Tint,?);',ic);   c_eval('Pi?.units = ''nPa'';', ic);
    c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',Tint, ?);',ic);
    % 2.3. fields
    c_eval('gseB?=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);
    % 2.4. Load defatt files
    c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',ic);
    c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',ic);

%%  2. load fast data, choose brst or fast
    % 2.1. PSD
    c_eval('efile? = mms.get_filepath([''mms?_fpi_fast_l2_des-dist''], time);', ic);
    c_eval('[ePDist?, ePDistError?] = mms.make_pdist(efile?);', ic); 
    
    
    c_eval('ifile? = mms.get_filepath([''mms?_fpi_fast_l2_dis-dist''], time);', ic);
    c_eval('[iPDist?, iPDistError?] = mms.make_pdist(ifile?);', ic); 
    % 2.2. moments
    c_eval('Tint = irf.tint(iPDist?.time.start, iPDist?.time.stop);', ic); Tint = Tint + [10 -10];
    c_eval('ni? = mms.get_data(''Ni_fpi_fast_l2'',Tint, ?);',ic);
    c_eval('Vi? = mms.get_data(''Vi_dbcs_fpi_fast_l2'',Tint,?);',ic);
    c_eval('Pi? = mms.get_data(''Pi_dbcs_fpi_fast_l2'',Tint,?);',ic);   c_eval('Pi?.units = ''nPa'';', ic);
    c_eval('ne? = mms.get_data(''Ne_fpi_fast_l2'',Tint, ?);',ic);
    
    %c_eval('ifile? = mms.get_filepath([''mms?_fpi_fast_l2_dis-moms''], time);', ic);
    %c_eval('tmpDataObj? = dataobj(ifile?);',ic);
    %c_eval(tmp,ic);
    %c_eval('ne? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_numberdensity_brst''));',ic);


    % 2.3. fields    
    c_eval('gseB?=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);
    c_eval('gseE?=mms.get_data(''E_gse_edp_fast_l2'', Tint, ?);',ic);
    % 2.4. Load defatt files
    c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',ic);
    c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',ic);

%%  3. compute
    % compute deHoffmanTeller frame
    [vht,eht,dvht,p,cc]=irf_vht(gseE1.tlim(tint_walen),gseB1.tlim(tint_walen));

    % 3.1. rotate to GSE
    c_eval('Pi? = mms.rotate_tensor(Pi?,''gse'',?);',ic);    
    c_eval('Vi? = mms_dsl2gse(Vi?,defatt?);',ic);
    c_eval('Pifac? = mms.rotate_tensor(Pi?,''fac'', gseB?);',ic);
 
    % 3.2. alpha: pressure anisotropy factor
    c_eval('alpha = irf_pres_anis(Pifac?, gseB?);', ic);
    
    % 3.3. gse to new123
    c_eval('B123 = irf_newxyz(gseB?, vec1, vec2, vec3);', ic);
    c_eval('Vi123 = irf_newxyz(Vi?, vec1, vec2, vec3);', ic);
    
    % 3.4. reference (MSH) region; in New frame (123);
    Bref = B123.tlim(tint_ref);                     Bref = irf.nanmean(Bref.data, 1);
    Viref = Vi123.tlim(tint_ref);                   Viref = irf.nanmean(Viref.data, 1);
    c_eval('Niref = ni?.tlim(tint_ref);', ic);      Niref = irf.nanmean(Niref.data, 1);
    alpharef = alpha.tlim(tint_ref);                alpharef = irf.nanmean(alpharef.data, 1);
    
    % 3.5. Vipred1: delta_B / sqrt(rho1)
    c_eval('B123 = B123.resample(ni?);', ic);
    c_eval('Vi123 = Vi123.resample(ni?);', ic);
    tmp1 = (B123 - Bref) * 21.8 / sqrt(Niref);
    Vipred1 = tmp1.resample(Vi123) * Jsign + Viref;
    
    % 3.6. Vipred2: B_2 / sqrt(rho2) - B_1 / sqrt(rho1)
    %c_eval('tmp2 = 21.8 * B123.data ./ sqrt([ni?.data ni?.data ni?.data]);', ic);
    %tmp2 = irf.ts_vec_xyz(B123.time, tmp2);
    %Vipred2 = (tmp2 - 21.8 * Bref / sqrt(Niref)) * Jsign + Viref;
    c_eval('tmp2 = 21.8 * (1-alpha) * B123 / sqrt(Niref) / sqrt(1-alpharef);', ic);       % Phan et al. [AG, 2004]
    Vipred2 = (tmp2 - 21.8 * sqrt(1-alpharef) * Bref / sqrt(Niref)) * Jsign + Viref;
    
    % 3.7. Vipred2: sqrt(1-alpha_2) * B_2 / sqrt(rho2) - sqrt(1-alpha_2) * B_1 / sqrt(rho1)  
    c_eval('tmp3 = 21.8 * sqrt(1-alpha) * B123 / sqrt(ni?);', ic);
    Vipred3 = (tmp3 - sqrt(1 - alpharef) * 21.8 * Bref / sqrt(Niref)) * Jsign + Viref;

    % 3.8. slope & CC
    Vi123w = Vi123.tlim(tint_walen);             Vipredw1 = Vipred1.tlim(tint_walen);    
    Vipredw2 = Vipred2.tlim(tint_walen);    Vipredw3 = Vipred3.tlim(tint_walen);    
    c_eval('p? = polyfit(Vipredw2.?.data, Vi123w.?.data, 1);     slope2? = p?(1);', ['x', 'y', 'z']);
    c_eval('corr? = corrcoef(Vipredw2.?.data, Vi123w.?.data);    cc2? = corr?(1, 2);', ['x', 'y', 'z']);    
   
%%  4. plot
    % 4.0. figure setting
    h = irf_plot(8,'newfigure');
    xSize=800; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.129;
    %c_eval('set(h(?), ''position'', [0.12 0.960-? * ywidth xwidth ywidth]);', 1: 7);
    
    isub = 1;
    % 4.1. Bxyz
    hca = irf_panel('Bgse');
    c_eval('irf_plot(hca, gseB?, ''LineWidth'', 1.0);', ic);
    ylabel(hca,{'B','[nT]'},'Interpreter','tex');
    irf_legend(hca,{'B_{x} ',' B_{y} ',' B_{z}'},[0.88 0.15])
            
    % 4.1. Vi
    hca = irf_panel('Vgse');
    c_eval('irf_plot(hca, Vi?, ''LineWidth'', 1.0);', ic);
    ylabel(hca,{'V_i','[km/s]'},'Interpreter','tex');
    irf_legend(hca,{'V_{i,x} ',' V_{i,y} ',' V_{i,z}'},[0.88 0.15])   
    
    % 4.2. Ni & Ne
    hca = irf_panel('ni');
    c_eval('irf_plot(hca, ni?, ''LineWidth'', 1.0);', ic);
    hold(hca, 'on');
    c_eval('irf_plot(hca, ne?, ''LineWidth'', 1.0);', ic);
    hold(hca, 'off');
    ylabel(hca,{'N','[cm^{-3}]'},'Interpreter','tex');
    irf_legend(hca,{'N_{i} ',' N_{e} '},[0.88 0.15])    
        
    % 4.3. ion energy flux
    hca = irf_panel('iEnflux');
    c_eval('irf_spectrogram(hca, iPDist?.deflux.omni.specrec, ''log'');', ic);
    hca.YScale = 'log';
    hca.YTick = 10.^[1 2 3 4];
    irf_legend(hca,'(c)',[0.99 0.98],'color','k')
    ylabel(hca,{'E_i','[eV]'},'Interpreter','tex');

    % 4.4. B123
    hca = irf_panel('B123');
    irf_plot(hca, B123, 'LineWidth', 1.0);
    ylabel(hca,{'B','[nT]'},'Interpreter','tex');
    irf_legend(hca,{'B_{1} ',' B_{2} ',' B_{3}'},[0.88 0.15])
    irf_legend(hca,'(d)',[0.99 0.98],'color','k')
    % vectors legend
    vec1_str = ['[', num2str(vec1(1),'% .2f'), ',', num2str(vec1(2),'% .2f'), ',', num2str(vec1(3),'% .2f'), ']'];
    vec2_str = ['[', num2str(vec2(1),'% .2f'), ',', num2str(vec2(2),'% .2f'), ',', num2str(vec2(3),'% .2f'), ']'];
    vec3_str = ['[', num2str(vec3(1),'% .2f'), ',', num2str(vec3(2),'% .2f'), ',', num2str(vec3(3),'% .2f'), ']'];
    irf_legend(hca, vec1_str, [1.003 0.8],'color','k', 'FontSize', 11)
    irf_legend(hca, vec2_str, [1.003 0.5],'color','b', 'FontSize', 11)
    irf_legend(hca, vec3_str, [1.003 0.1],'color','r', 'FontSize', 11)
   
    % 4.5. Vi123.x
    hca = irf_panel('Vi1');
    irf_plot(hca, Vipredw2.x, 'LineWidth', 1.0, 'color', 'r');
    hold(hca, 'on');
    irf_plot(hca, Vi123.x, 'LineWidth', 1.0, 'color', 'k');
    irf_plot(hca, Vipred2.x, '--', 'LineWidth', 1.0, 'color', 'r');    
    hold(hca, 'off');
    ylabel(hca,{'V_1','[km/s]'},'Interpreter','tex');
    irf_legend(hca,{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster');
    irf_legend(hca,'(e)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2x, '% .2f')];
    ccx_str = ['cc=', num2str(cc2x, '% .2f')];
    irf_legend(hca, slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(hca, ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.6. Vi123.y
    hca = irf_panel('Vi2');
    irf_plot(hca, Vipredw2.y, 'LineWidth', 1.0, 'color', 'r');
    hold(hca, 'on');
    irf_plot(hca, Vi123.y, 'LineWidth', 1.0, 'color', 'k');
    irf_plot(hca, Vipred2.y, '--', 'LineWidth', 1.0, 'color', 'r');    
    hold(hca, 'off');
    ylabel(hca,{'V_2','[km/s]'},'Interpreter','tex');
    irf_legend(hca,{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster')
    irf_legend(hca,'(f)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2y, '% .2f')];
    ccx_str = ['cc=', num2str(cc2y, '% .2f')];
    irf_legend(hca, slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(hca, ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.5. Vi123.z
    hca = irf_panel('Vi3');
    irf_plot(hca, Vipredw2.z, 'LineWidth', 1.0, 'color', 'r');
    hold(hca, 'on');
    irf_plot(hca, Vi123.z, 'LineWidth', 1.0, 'color', 'k');
    irf_plot(hca, Vipred2.z, '--', 'LineWidth', 1.0, 'color', 'r');    
    hold(hca, 'off');
    ylabel(hca,{'V_3','[km/s]'},'Interpreter','tex');
    irf_legend(hca,{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster')
    irf_legend(hca,'(g)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2z, '% .2f')];
    ccx_str = ['cc=', num2str(cc2z, '% .2f')];
    irf_legend(hca, slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(hca, ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.X. global
    load('caa/cmap.mat');
    colormap(irf_panel('iEnflux'),cmap); 
    title(h(1), strcat('MMS ', num2str(ic)));
    irf_plot_axis_align(h)
    irf_zoom(h(6:8),'y')
    irf_zoom(h,'x', tint);
    irf_pl_mark(h(6:8), tint_ref, 'red')
    irf_pl_mark(h(6:8), tint_walen, 'yellow'); %, 'LineStyle', '--', 'LineWidth', 1.0)    

%%

