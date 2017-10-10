%% Load data
ic = 1:4;
tint = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z'); %20151112071854


c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst'',tint,?);',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4, c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else, c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end
%% Prepare data
units = irf_units; 
gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km';
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';
avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
gseOhmGradPe = gseGradPe/avNe.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseOhmGradPe.units = 'mV/m';

% Currents from moments
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 
c_eval('[gseJ?par,gseJ?perp] = irf_dec_parperp(gseB?,gseJ?); gseJ?par.name = ''J par''; gseJ?perp.name = ''J perp'';',ic)

c_eval('dbcsJe? = -units.e*ne?*dbcsVe?*1e3*1e6*1e9; dbcsJe?.units = ''nA/m^2''; dbcsJe?.coordinateSystem = ''DBCS'';',ic);
c_eval('dbcsJi? = units.e*ne?*dbcsVi?.resample(ne?.time)*1e3*1e6*1e9; dbcsJi?.units = ''nA/m^2''; dbcsJi?.coordinateSystem = ''DBCS'';',ic);
c_eval('dbcsJ? = (dbcsJe?+dbcsJi?);',ic);
dbcsAvJ = (dbcsJ1+dbcsJ2.resample(dbcsJ1.time)+dbcsJ3.resample(dbcsJ1.time)+dbcsJ4.resample(dbcsJ1.time))/4; 
c_eval('[dbcsJ?par,dbcsJ?perp] = irf_dec_parperp(dmpaB?,dbcsJ?); dbcsJ?par.name = ''J par''; dbcsJ?perp.name = ''J perp'';',ic)

c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)
c_eval('dbcsVexB? = dbcsVe?.cross(dmpaB?.resample(dbcsVe?))*1e-3; dbcsVexB?.units = ''mV/m'';',ic)
c_eval('dbcsVixB? = dbcsVi?.cross(dmpaB?.resample(dbcsVi?))*1e-3; dbcsVixB?.units = ''mV/m'';',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[dbcsVe?par,dbcsVe?perp] = irf_dec_parperp(dmpaB?,dbcsVe?); dbcsVe?par.name = ''Ve par''; dbcsVe?perp.name = ''Ve perp'';',ic)

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(dmpaB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('[dslE?par,dslE?perp] = irf_dec_parperp(dmpaB?,dslE?); dslE?par.name = ''E par''; dslE?perp.name = ''E perp'';',ic)

%% Reconstruct Electric field from probe voltages
baselength = 2*14.5;
multiplier = 1.5;
offsetE1z = -2.5;
offsetE2z = -3.5;
offsetE3z = 0;
offsetE4z = 0.5;
c_eval('dslE?z = irf.ts_scalar(dcv?.time,(dcv?.data(:,6)-dcv?.data(:,5))/baselength*multiplier*1e3-offsetE?z);',1:4)
c_eval('dslE?_newz = dslE?.clone(dslE?.time,[dslE?.data(:,1:2) dslE?z.data]);',1:4)
%%
% lowpass filter E field and subtract this part
ffilt = 0.15;
c_eval('[dslE?par_newz,dslE?perp_newz] = irf_dec_parperp(dmpaB?,dslE?_newz); dslE?par_newz.units =''mv/m''; dslE?perp_newz.units =''mv/m'';',ic)
%c_eval('dslE?_lowf = dslE?perp_newz.filt(0,ffilt,[],3);',ic)
c_eval('dslE?_lowf = dslE?_newz.filt(0,ffilt,[],3);',ic)
c_eval('dbcsVexB?_lowf = dbcsVexB?.filt(0,ffilt,[],3);',ic)
c_eval('E?_vefit = dslE?_lowf + dbcsVexB?_lowf.resample(dslE?);',ic)
c_eval('dslE?_detrend = dslE?_newz - E?_vefit;',ic);
c_eval('dslE?_new = dslE?_detrend; dslE?_new.name = ''E new'';',ic)

c_eval('[dslE?par_new,dslE?perp_new] = irf_dec_parperp(dmpaB?,dslE?_new); dslE?par_new.name = ''E par''; dslE?perp_new.name = ''E perp'';',ic)
%% Transform new E field to gse coordinte system
% c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',1:4)
% c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',1:4)
load /Users/Cecilia/Data/MMS/2015Nov12/defatt.mat
c_eval('gseE?_new = mms_dsl2gse(dslE?_new,defatt?);',ic)
c_eval('[gseE?par_new,gseE?perp_new] = irf_dec_parperp(gseB?,gseE?_new); gseE?par_new.name = ''E par''; gseE?perp_new.name = ''E perp'';',ic)

c_eval('gseEVexB?_new = gseE?_new.resample(gseVexB?.time) + gseVexB?; gseEVexB?_new.name = ''E+VexB'';',ic)
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

gseAvE_new = 0.25*(gseE1_new + gseE2_new.resample(gseE1_new) + gseE3_new.resample(gseE1_new) + gseE4_new.resample(gseE1_new)); gseAvE_new.name = 'av E new';
gseAvVexB = 0.25*(gseVexB1 + gseVexB2.resample(gseVexB1) + gseVexB3.resample(gseVexB1) + gseVexB4.resample(gseVexB1)); gseVexB.name = 'av Ve new';

%% Transform new E field to LMN (mva) coordinte system
tint_bcs_utc = '2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z';
tint_bcs = irf.tint(tint_bcs_utc);
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(tint_bcs));',1:4)
c_eval('mva_all(:,:,?) = v?;',1:4)
mva_mean = mean(mva_all,3);L = mva_mean(1,:); M = mva_mean(2,:); N = mva_mean(3,:);
coordLabels = {'L','M','N'};
lmn = [L;M;N];
c_eval('mvaE?_new = gseE?_new*lmn'';',ic)
c_eval('mvaE?par_new = gseE?par_new; mvaE?perp_new = gseE?perp_new*lmn'';',ic)
mvaAvE_new = gseAvE_new*lmn'';
c_eval('mvaEVexB?_new = gseEVexB?_new*lmn'';',ic)

%% E dot J, E+VexB fot J etc

 %c_eval('dbcsVe?par.tlim(tint)*(dslE?par.resample(dbcsVe?par));',ic)

c_eval('EdotJ? = gseE?.resample(gseJ?).dot(gseJ?)/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?_new = gseE?_new.resample(gseJ?).dot(gseJ?)/1000; EdotJ?_new.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?par_new = gseE?par_new.resample(gseJ?par)*gseJ?par/1000; EdotJ?par_new.units = ''nW/m^3''; EdotJ?par_new.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?perp_new = gseE?perp_new.resample(gseJ?perp).dot(gseJ?perp)/1000; EdotJ?perp_new.units = ''nW/m^3''; EdotJ?perp_new.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?par = gseE?par.resample(gseJ?par)*gseJ?par/1000; EdotJ?par.units = ''nW/m^3''; EdotJ?par.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?perp = gseE?perp.resample(gseJ?perp).dot(gseJ?perp)/1000; EdotJ?perp.units = ''nW/m^3''; EdotJ?perp.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?parL3 = E?par.y.resample(gseJ?par)*gseJ?par/1000; EdotJ?parL3.units = ''nW/m^3''; EdotJ?parL3.name = ''E*J par'';',[1 4]); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('RedotJ? = gseEVexB?.resample(gseJ?).dot(gseJ?)/1000; RedotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('RedotJ?_new = gseEVexB?_new.resample(gseJ?).dot(gseJ?)/1000; RedotJ?_new.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

%% Figure: 1 sc overview electric field
tintZoom = irf.tint('2015-11-12T07:19:18.00Z/2015-11-12T07:19:28.00Z');
npanels = 9;
cmap = 'jet';

% Initialize figure
ic = 4;
h = irf_plot(npanels);

% Plot
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 1 % Vi dsl/dbcs/dmpa X
  hca = irf_panel('E dsl X');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{dslE?.tlim(tint).x,-dbcsVexB?.x.tlim(tint),-dbcsVixB?.x.tlim(tint),dslE?_newz.tlim(tint).x,gseE?.tlim(tint).x},''comp'');',ic)
  c_eval('irf_plot(hca,{dslE?.tlim(tint).x,-dbcsVexB?.x.tlim(tint),-dbcsVixB?.x.tlim(tint),irf.ts_vec_xyz(dslE?.time,dslE?.data*NaN).x.tlim(tint),dslE?_detrend.tlim(tint).x},''comp'');',ic)
  hca.YLabel.String = {'E_x','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'E','-v_exB','-v_ixB','E_{new z}','E_{GSE}'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'E','-v_exB','-v_ixB','','E_{detrend}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Vi dsl/dbcs/dmpa Y
  hca = irf_panel('E dsl Y');
  set(hca,'ColorOrder',mms_colors('matlab'))  
  %c_eval('irf_plot(hca,{dslE?.tlim(tint).y,-dbcsVexB?.y.tlim(tint),-dbcsVixB?.y.tlim(tint),dslE?_newz.tlim(tint).y,gseE?.tlim(tint).y},''comp'');',ic)
  c_eval('irf_plot(hca,{dslE?.tlim(tint).y,-dbcsVexB?.y.tlim(tint),-dbcsVixB?.y.tlim(tint),irf.ts_vec_xyz(dslE?.time,dslE?.data*NaN).y.tlim(tint),dslE?_detrend.tlim(tint).y},''comp'');',ic)
  hca.YLabel.String = {'E_y','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'E','-v_exB','-v_ixB','E_{new z}','E_{GSE}'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'E','-v_exB','-v_ixB','','E_{detrend}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Vi dsl/dbcs/dmpa X
  hca = irf_panel('E dsl z');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{dslE?.tlim(tint).z,-dbcsVexB?.z.tlim(tint),-dbcsVixB?.z.tlim(tint),dslE?_newz.tlim(tint).z,gseE?.tlim(tint).z},''comp'');',ic)
  c_eval('irf_plot(hca,{dslE?.tlim(tint).z,-dbcsVexB?.z.tlim(tint),-dbcsVixB?.z.tlim(tint),dslE?_newz.tlim(tint).z,dslE?_detrend.tlim(tint).z},''comp'');',ic)
  hca.YLabel.String = {'E_z','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'E','-v_exB','-v_ixB','E_{new z}','E_{GSE}'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'E','-v_exB','-v_ixB','E_{new z}','E_{detrend}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Electric field subtracted to detrend
  hca = irf_panel('E dsl detrend');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{E?_vefit.tlim(tint).x,E?_vefit.tlim(tint).y,E?_vefit.tlim(tint).z},''comp'');',ic)
  hca.YLabel.String = {'E+v_exB offs','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{sprintf('f_{lowpass} = %g Hz',ffilt)},[0.05 0.9],'fontsize',12,'color',[0 0 0]);  
end
if 1 % Electric field + VexB new
  hca = irf_panel('E + VexB new');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dslE?perp_new.tlim(tint).x+dbcsVexB?.resample(dslE?_new).x,dslE?perp_new.tlim(tint).y+dbcsVexB?.resample(dslE?_new).y,dslE?perp_new.tlim(tint).z+dbcsVexB?.resample(dslE?_newz).z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}+v_exB','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'new'},[0.05 0.9],'fontsize',12,'color',[0 0 0]);  
end
if 1 % Electric field + VexB old
  hca = irf_panel('E + VexB old');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dslE?perp.tlim(tint).x+dbcsVexB?.resample(dslE?).x,dslE?perp.tlim(tint).y+dbcsVexB?.resample(dslE?).y,dslE?perp.tlim(tint).z+dbcsVexB?.resample(dslE?).z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}+v_exB','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12); 
  irf_legend(hca,{'original'},[0.05 0.9],'fontsize',12,'color',[0 0 0]);  
end
if 1 % E par new and old
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dslE?par.tlim(tint),dslE?par_new.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);   
end
if 1 % gradPe Ohm (mV/m)
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{-gseOhmGradPe.x,-gseOhmGradPe.y,-gseOhmGradPe.z},'comp');
  hca.YLabel.String = {'-\nabla \cdot P_e','(GSE)','(mV/m)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'X','Y','Z'},[0.98 0.9],'fontsize',12);
end
irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).YLim = [-4 6];',6:7);
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure: 1 sc overview electric field, incl E par and E dot J
tintZoom = irf.tint('2015-11-12T07:19:18.00Z/2015-11-12T07:19:28.00Z');
npanels = 9;
cmap = 'jet';

% Initialize figure
ic = 4;
h = irf_plot(npanels);

% Plot
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E new
  hca = irf_panel('E new');
  set(hca,'ColorOrder',mms_colors('xyz'))
 
  c_eval('irf_plot(hca,{dslE?_new.x.tlim(tint),dslE?_new.y.tlim(tint),dslE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % J new
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyz')) 
  c_eval('irf_plot(hca,{dbcsJ?.x.tlim(tint),dbcsJ?.y.tlim(tint),dbcsJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(DBCS)','(nT/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E new
  hca = irf_panel('E new perp par');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  c_eval('irf_plot(hca,{dslE?perp_new.x.tlim(tint),dslE?perp_new.y.tlim(tint),dslE?perp_new.z.tlim(tint),dslE?par_new.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);  
end
if 1 % Electric field + VexB new
  hca = irf_panel('E + VexB new');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  %c_eval('irf_plot(hca,{dslE?perp_new.tlim(tint).x+dbcsVexB?.resample(dslE?_new).x,dslE?perp_new.tlim(tint).y+dbcsVexB?.resample(dslE?_new).y,dslE?perp_new.tlim(tint).z+dbcsVexB?.resample(dslE?_newz).z},''comp'');',ic)
  c_eval('irf_plot(hca,{dslE?_new.tlim(tint).x+dbcsVexB?.resample(dslE?_new).x,dslE?_new.tlim(tint).y+dbcsVexB?.resample(dslE?_new).y,dslE?_new.tlim(tint).z+dbcsVexB?.resample(dslE?_newz).z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}+v_exB','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'new'},[0.05 0.9],'fontsize',12,'color',[0 0 0]);  
end
if 0 % Electric field + VexB old
  hca = irf_panel('E + VexB old');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dslE?perp.tlim(tint).x+dbcsVexB?.resample(dslE?).x,dslE?perp.tlim(tint).y+dbcsVexB?.resample(dslE?).y,dslE?perp.tlim(tint).z+dbcsVexB?.resample(dslE?).z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}+v_exB','(DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12); 
  irf_legend(hca,{'original'},[0.05 0.9],'fontsize',12,'color',[0 0 0]);  
end
if 1 % E par new and old
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  try
    c_eval('irf_plot(hca,{dslE?par.tlim(tint),dslE?par_new.tlim(tint),E?par.y},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('xyz'))  
    irf_legend(hca,{'original','new','L3'},[0.98 0.9],'fontsize',12);   
  catch
    c_eval('irf_plot(hca,{dslE?par.tlim(tint),dslE?par_new.tlim(tint)},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('xyz'))  
    irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);   
  end
  hca.YLabel.String = {'E_{||}','(DSL)','(mV/m)'};  
end
if 1 % J par perp 
  hca = irf_panel('J perp par');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  c_eval('irf_plot(hca,{dbcsJ?perp.x.tlim(tint),dbcsJ?perp.y.tlim(tint),dbcsJ?perp.z.tlim(tint),dbcsJ?par.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(DBCS)','(nT/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);  
end
if 0 % J par
  hca = irf_panel('J par');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dbcsJ?par.tlim(tint)},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyz'))  
  hca.YLabel.String = {'J_{||}','(DBCS)','(nA/m^2)'};  
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{dbcsVe?par.tlim(tint)},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyz'))  
  hca.YLabel.String = {'V_{e,||}','(DBCS)','(mV/m)'};  
end
if 1 % E dot J
  hca = irf_panel('E dot J');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{EdotJ?.tlim(tint),EdotJ?_new},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);      
  hca.YLabel.String = {'E \cdot J','(DSL/DBCS)','(nW/m^3)'};  
end
if 1 % E dot J par
  hca = irf_panel('E dot J par');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  try    
    c_eval('irf_plot(hca,{EdotJ?par,EdotJ?par_new.tlim(tint),EdotJ?parL3.tlim(tint)},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('xyz'))  
    irf_legend(hca,{'original','new','L3'},[0.98 0.9],'fontsize',12);  
  catch
    c_eval('irf_plot(hca,{EdotJ?par,EdotJ?par_new.tlim(tint)},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('xyz'))  
    irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);  
  end  
  set(hca,'ColorOrder',mms_colors('xyz'))  
  hca.YLabel.String = {'E_{||} \cdot J_{||}','(DSL/DBCS)','(nW/m^3)'};  
end
if 1 % E dot J perp
  hca = irf_panel('E dot J perp');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{EdotJ?perp.tlim(tint),EdotJ?perp_new},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);      
  hca.YLabel.String = {'E_{\perp} \cdot J_{\perp}','(DSL/DBCS)','(nW/m^2)'};  
end
if 1 % (E + VexB) dot J
  hca = irf_panel('(E+VexB) dot J');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  c_eval('irf_plot(hca,{RedotJ?.tlim(tint),RedotJ?_new},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'original','new'},[0.98 0.9],'fontsize',12);      
  hca.YLabel.String = {'(E+v_exB) \cdot J','(DSL/DBCS)','(nW/m^3)'};  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%c_eval('h(?).YLim = [-4 6];',6:7);
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end


%% Figure 1: 1 sc overview electric field
%tint = irf.tint('2015-10-16T10:34:20.00Z/2015-10-16T10:35:30.00Z');
npanels = 11;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
if 1 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 0 % B mva
  hca = irf_panel('B mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 0 % Vi mva
  hca = irf_panel('Vi mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint),mvaVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % Vi gse
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.deflux.tlim(tint).pitchangles(dmpaB?,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 0 % Ve mva
  hca = irf_panel('Ve mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','M'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve perp par mva
  hca = irf_panel('Ve perp par mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L,\perp','M,\perp','N,\perp','||'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve gse
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve perp par gse
  hca = irf_panel('Ve perp par gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x,\perp','y,\perp','z,\perp','||'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % V ExB LMN
  hca = irf_panel('V ExB mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x.tlim(tint),mvaVExB?.y.tlim(tint),mvaVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 1 % V ExB GSE
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.tlim(tint),gseVExB?.y.tlim(tint),gseVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
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
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
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
  colormap(hca,cmap)   
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 200];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa 64
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 1000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % E mva
  hca = irf_panel('E mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E par
  hca = irf_panel('E par');
  try
    set(hca,'ColorOrder',mms_colors('b12'))
    c_eval('irf_plot(hca,{E?par.x,E?par.y,gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('b12'))
    irf_legend(hca,{'error','E_{||}','E_{||}'},[0.98 0.9],'fontsize',12);
  catch
    set(hca,'ColorOrder',mms_colors('2'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);
  end
  %irf_zoom(hca,'y')
end
if 0 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.x.tlim(tint),mvaE2.x.tlim(tint),mvaE3.x.tlim(tint),mvaE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{L}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.y.tlim(tint),mvaE2.y.tlim(tint),mvaE3.y.tlim(tint),mvaE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{M}','(mV/m)'};
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.z.tlim(tint),mvaE2.z.tlim(tint),mvaE3.z.tlim(tint),mvaE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{N}','(mV/m)'};
end
if 1 % E DSL X
  hca = irf_panel('E DSL X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{X,DSL}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % E DSL Y
  hca = irf_panel('E DSL Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.y.tlim(tint),dslE2.y.tlim(tint),dslE3.y.tlim(tint),dslE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Y,DSL}','(mV/m)'};
end
if 1 % E DSL Z
  hca = irf_panel('E DSL Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.z.tlim(tint),dslE2.z.tlim(tint),dslE3.z.tlim(tint),dslE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Z,DSL}','(mV/m)'};
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law N
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm LMN residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = mvaAvE + mvaOhmVexB.resample(mvaAvE) + mvaOhmGradPe.resample(mvaAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
end
if 1 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm X');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  irf_plot(hca,{gseAvE.x,gseOhmGradPe.x,-gseOhmVexB.x},'comp');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E_Z','(mV/m)'};
end
if 1 % Ohm's law Y
  hca = irf_panel('Ohm Y');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.y,gseOhmGradPe.y,-gseOhmVexB.y},'comp');
  hca.YLabel.String = {'E_Y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law Z
  hca = irf_panel('Ohm Z');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.z,gseOhmGradPe.z,-gseOhmVexB.z},'comp');
  hca.YLabel.String = {'E_Z','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm XYZ residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = gseAvE + gseOhmVexB.resample(gseAvE) + gseOhmGradPe.resample(gseAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'X','Y','Z'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%% Figure 1: 1 sc overview electric field, spacecraft coord. sys.
%tint = irf.tint('2015-10-16T10:34:20.00Z/2015-10-16T10:35:30.00Z');
npanels = 9;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B_{DMPA}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{scPot1*(-1),scPot2*(-1),scPot3*(-1),scPot4*(-1)},'comp');
  hca.YLabel.String = {'-V_{SC}','(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))    
end
if 0 % Vi gse
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % Ve perp par mva
  hca = irf_panel('Ve perp par mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L,\perp','M,\perp','N,\perp','||'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve gse
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Ve perp par gse
  hca = irf_panel('Ve perp par gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x,\perp','y,\perp','z,\perp','||'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % V ExB
  hca = irf_panel('V ExB mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x.tlim(tint),mvaVExB?.y.tlim(tint),mvaVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 0 % V ExB
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.tlim(tint),gseVExB?.y.tlim(tint),gseVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
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
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % E par
  hca = irf_panel('E par');
  try
    set(hca,'ColorOrder',mms_colors('b12'))
    c_eval('irf_plot(hca,{E?par.x,E?par.y,gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('b12'))
    irf_legend(hca,{'error','E_{||}','E_{||}'},[0.98 0.9],'fontsize',12);
  catch
    set(hca,'ColorOrder',mms_colors('2'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);
  end
  %irf_zoom(hca,'y')
end
if 1 % E DSL X
  hca = irf_panel('E DSL X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{X,DSL}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % E DSL Y
  hca = irf_panel('E DSL Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.y.tlim(tint),dslE2.y.tlim(tint),dslE3.y.tlim(tint),dslE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Y,DSL}','(mV/m)'};
end
if 1 % E DSL Z
  hca = irf_panel('E DSL Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.z.tlim(tint),dslE2.z.tlim(tint),dslE3.z.tlim(tint),dslE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Z,DSL}','(mV/m)'};
end
if 1 % ve x B GSE X
  hca = irf_panel('Ve x B GSE X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.x.tlim(tint),-gseVexB2.x.tlim(tint),-gseVexB3.x.tlim(tint),-gseVexB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{X,GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve x B GSE Y
  hca = irf_panel('Ve x B GSE Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.y.tlim(tint),-gseVexB2.y.tlim(tint),-gseVexB3.y.tlim(tint),-gseVexB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{Y,GSE}','(mV/m)'};
end
if 1 % Ve x B GSE Z
  hca = irf_panel('Ve x B GSE Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.z.tlim(tint),-gseVexB2.z.tlim(tint),-gseVexB3.z.tlim(tint),-gseVexB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{Z,GSE}','(mV/m)'};
end
if 1 % gradPe GSE
  hca = irf_panel('gradPe/ne gse mv/m');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseOhmGradPe.x,gseOhmGradPe.y,gseOhmGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e/ne','(mv/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % Ohm's law x
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law y
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y,-mvaOhmVixB.y,mvaOhmJxB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law z
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z,-mvaOhmVixB.z,mvaOhmJxB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm LMN residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = mvaAvE + mvaOhmVexB.resample(mvaAvE) + mvaOhmGradPe.resample(mvaAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
end

c_eval('h(?).YLim = [-4 4];',3:9)
irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%% Figure 2: 4 sc overview electric field potentials, spacecraft coord. sys.
tint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:35.00Z');
npanels = 10;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
ic = 1;
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,irf_ssub('MMS ?',ic),[0.02,0.1])
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{scPot1*(-1),scPot2*(-1),scPot3*(-1),scPot4*(-1)},'comp');  
  hca.YLabel.String = {irf_ssub('-V_{SC}',ic),'(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))    
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.1],'fontsize',12);
end
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 2;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 3;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 4;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end


irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure 2: 1 sc overview electric field potentials, spacecraft coord. sys.
%tint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:35.00Z');
npanels = 5;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
ic = 4;
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,scPot1*(-1));',ic)  
  ylabel(hca,{irf_ssub('-V_{SC}',ic),'(V)'},'interpreter','tex')
  %hca.YLabel.String = {irf_ssub('-V_{SC}',ic),'(V)'};    
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  %set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  ylabel(hca,{'DCV','(V)'},'interpreter','tex')
  irf_legend(hca,{'1','2','3','4','5','6'},[0.98 0.9],'fontsize',12);  
  %set(hca,'ColorOrder',mms_colors('1234'))    
end
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? (DSL)',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % difference between axial probes
  hca = irf_panel(irf_ssub('V6-V5',ic));
  baselength = 14.5;
  c_eval('axdiff = irf.ts_scalar(dcv?.time,(dcv?.data(:,6)-dcv?.data(:,5))/baselength*1e3);',ic)
  set(hca,'ColorOrder',mms_colors('z'))
  c_eval('irf_plot(hca,axdiff);',ic)
  ylabel(hca,{sprintf('(V_6 - V_5)/%.1f',baselength),'(mV/m)'},'interpreter','tex')
  %hca.YLabel.String = {sprintf('(V_6 - V_5)/%.1f',baselength),'(mV/m)'};    
end

h(1).Title.String = irf_ssub('MMS ?',ic);

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end
