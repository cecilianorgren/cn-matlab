ic = 1:4;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tintzoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
localuser = datastore('local','user');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

%% Load data
ic = 1:4;
units = irf_units;
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('E?par = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% EDI
c_eval('jedi?  = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)
c_eval('jedi?_err = mms.get_data(''Flux-err-amb-pm2_edi_brst_l2'',tint,?);',ic)
c_eval(['jedi?_err_plus   = jedi?+jedi?_err;' ...  
        'jedi?_err_minus  = jedi?+(jedi?_err*-1);'],ic)

if 0
  %%
  c_eval('jedi? = ePitch?_flux_edi;',ic)
  c_eval('jedi?_err = ePitch?_flux_edi_err;',ic)
  c_eval('jedi?_err_minus = ePitch?_flux_edi_err_minus;',ic)
  c_eval('jedi?_err_plus = ePitch?_flux_edi_err_plus;',ic)
  
end
      
% div E

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

%% Figure
npanels = 4;
h = irf_plot(npanels);

if 1 % Epar 1234
  hca = irf_panel('Epar 1234');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par.tlim(tintzoom),gseE2par.tlim(tintzoom),gseE3par.tlim(tintzoom),gseE4par.tlim(tintzoom)},'comp');
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % EDI 1234
  hca = irf_panel('EDI 1234');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1.palim(pa).tlim(tintzoom),jedi2.palim(pa).tlim(tintzoom),jedi3.palim(pa).tlim(tintzoom),jedi4.palim(pa).tlim(tintzoom)},'comp');
  hca.YLabel.String = {'j^{EDI}','(cm^{-2} s^{-1} sr^{-1})'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % EDI 1234
  hca = irf_panel('EDI err 1234');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1_err.palim(pa).tlim(tintzoom),jedi2_err.palim(pa).tlim(tintzoom),jedi3_err.palim(pa).tlim(tintzoom),jedi4_err.palim(pa).tlim(tintzoom)},'comp');
  hca.YLabel.String = {'j^{EDI} (err)','(cm^{-2} s^{-1} sr^{-1})'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 0 % EDI err 1
  hca = irf_panel('EDI err 1');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1.palim(pa),jedi1_err_minus.palim(pa),jedi1_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms1)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'j','j-j_{err}','j+j_{err}'},[0.98 0.98])
end
if 0 % EDI err 2
  hca = irf_panel('EDI err 2');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi2.palim(pa),jedi2_err_minus.palim(pa),jedi2_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms2)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'j','j-j_{err}','j+j_{err}'},[0.98 0.98])
end
if 0 % EDI err 3
  hca = irf_panel('EDI err 3');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi3.palim(pa),jedi3_err_minus.palim(pa),jedi3_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms3)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'j','j-j_{err}','j+j_{err}'},[0.98 0.98])
end
if 0 % EDI err 4
  hca = irf_panel('EDI err 4');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi4.palim(pa),jedi4_err_minus.palim(pa),jedi4_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms4)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'j','j-j_{err}','j+j_{err}'},[0.98 0.98])
end
if 1 % EDI err 1234, patch - pretty but quite cost intensive, so only do for tintzoom interval 
  hca = irf_panel('EDI err 1234 patch');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  hold(hca,'on')
  hp = irf_patch(hca,{jedi1_err_minus.palim(pa).tlim(tintzoom),jedi1_err_plus.palim(pa).tlim(tintzoom)});
  hp.FaceColor = mms_colors('1');
  hp.EdgeColor = mms_colors('1');
  hold(hca,'on')  
  hp = irf_patch(hca,{jedi2_err_minus.palim(pa).tlim(tintzoom),jedi2_err_plus.palim(pa).tlim(tintzoom)});
  hp.FaceColor = mms_colors('2');
  hp.EdgeColor = mms_colors('2');
  hold(hca,'on')
  hp = irf_patch(hca,{jedi3_err_minus.palim(pa).tlim(tintzoom),jedi3_err_plus.palim(pa).tlim(tintzoom)});
  hp.FaceColor = mms_colors('3');
  hp.EdgeColor = mms_colors('3');
  hold(hca,'on')
  hp = irf_patch(hca,{jedi4_err_minus.palim(pa).tlim(tintzoom),jedi4_err_plus.palim(pa).tlim(tintzoom)});
  hp.FaceColor = mms_colors('4');
  hp.EdgeColor = mms_colors('4');
  hold(hca,'off')
  %irf_plot(hca,{jedi1.palim(pa),jedi2.palim(pa).palim(pa),jedi3.palim(pa),jedi4.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI}','(cm^{-2} s^{-1} sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end

irf_zoom(h,'x',tintzoom)
%irf_zoom(h,'y')
%% Figure, timeseries and scatter
npanels = 3;
nrows = 3;
ncols = 2;
[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');

if 1 % EDI 1234
  hca = irf_panel('EDI 1234');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1.palim(pa),jedi2.palim(pa).palim(pa),jedi3.palim(pa),jedi4.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI}','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % EDI 1234
  hca = irf_panel('EDI err 1234');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1_err.palim(pa),jedi2_err.palim(pa),jedi3_err.palim(pa),jedi4_err.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (err)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 1 % EDI err 1
  hca = irf_panel('EDI err 1');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{jedi1.palim(pa),jedi1_err_minus.palim(pa),jedi1_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms1)','(cm^{-2} s^{-1})'};
  irf_legend(hca,{'j','j-j_{err}','j+j_{err}'},[0.98 0.98])
end
if 0 % EDI err 2
  hca = irf_panel('EDI err 2');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch2_flux_edi_err_minus.palim(pa),ePitch2_flux_edi_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms2)','(cm^{-2} s^{-1})'};
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 0 % EDI err 3
  hca = irf_panel('EDI err 3');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch3_flux_edi_err_minus.palim(pa),ePitch3_flux_edi_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms3)','(cm^{-2} s^{-1})'};
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end
if 0 % EDI err 4
  hca = irf_panel('EDI err 4');
  pa = [170 180];
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch4_flux_edi_err_minus.palim(pa),ePitch4_flux_edi_err_plus.palim(pa)},'comp');
  hca.YLabel.String = {'j^{EDI} (mms4)','(cm^{-2} s^{-1})'};
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.98])
end

irf_zoom(h,'x',tintzoom)
%%
nedges = 20;
jedges = linspace(0,1.5e7,nedges);
isub = 1;
if 1 % j1,j2
  hca = h2(isub); isub = isub + 1;
  A = jedi1.palim(pa).data;
  B = jedi2.palim(pa).resample(jedi1).data;
  [N,XEDGES,YEDGES] = histcounts2(A,B,jedges,jedges);
  surf(hca,jedges,jedges,zeros(nedges,nedges)',log10(N)');
  shading(hca,'flat')
  view(hca,[0 0 1])
end

