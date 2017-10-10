% See Scudder2012

% Make field aligned coordinates
tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
c_eval('[out?,l?,v?]=irf_minvar(dmpaB?.tlim(tint_mva));',1:4)
v = v3;

% Parallel (z): parallel to B/L
% Perp1 (x): close to M
% Perp2 (y): close to N
ic = 1:4;
c_eval('fsi = 1/(ni?brst.time(2)-ni?brst.time(1));',ic); fny = fsi/2;
c_eval('Pi?_lowres = irf_filt(Pi?brst,0,fny/2,fsi,5);',ic)
c_eval('Ti?_lowres = irf_filt(Ti?brst,0,fny/2,fsi,5);',ic)
c_eval('ni?_lowres = irf_filt(ni?brst,0,fny/2,fsi,5);',ic)

c_eval('fse = 1/(ne?brst.time(2)-ne?brst.time(1));',ic); fny = fse/2;
c_eval('Pe?_lowres = irf_filt(Pe?brst,0,fny/2,fse,5);',ic)
c_eval('Te?_lowres = irf_filt(Te?brst,0,fny/2,fse,5);',ic)
c_eval('ne?_lowres = irf_filt(ne?brst,0,fny/2,fse,5);',ic)



c_eval('par = dmpaB?/dmpaB?.abs;',ic);
c_eval('perp2 = par.cross(irf.ts_vec_xyz(par.time,repmat(v(2,:),par.length,1))); perp2 = perp2/perp2.abs;',ic);
c_eval('perp1 = perp2.cross(par); perp1 = perp1/perp1.abs;',ic);
par.units = ''; perp1.units = ''; perp2.units = '';

c_eval('facPe? = irf.ts_vec_xyz(Pe?_lowres.time,[Pe?_lowres.dot(perp1.resample(Pe?_lowres.time)).data,Pe?_lowres.dot(perp2.resample(Pe?_lowres.time)).data,Pe?_lowres.dot(par.resample(Pe?_lowres.time)).data]);',ic);
c_eval('facTe? = irf.ts_vec_xyz(Te?_lowres.time,[Te?_lowres.dot(perp1.resample(Te?_lowres.time)).data,Te?_lowres.dot(perp2.resample(Te?_lowres.time)).data,Te?_lowres.dot(par.resample(Te?_lowres.time)).data]);',ic);
%%
ic = 1;
% Define parameters
% Electron thermal anisotropy;
c_eval('Ane? = irf.ts_scalar(facTe?.time,facTe?.z.data./(facTe?.x.data+facTe?.y.data)/2);',ic);
% Agyrotropy of the measured electron pressure tensor
% Scudder: Agy = 2*|1-alpha|/(1-alpha) > 1; alpha = Pperp1/Pperp2;
% I think the following is better:
c_eval('Agy? = irf.ts_scalar(facPe?.time,log10(facPe?.x.data./facPe?.y.data));',ic);





%% Make plots
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
ic = 1;

h = irf_plot(9);

hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.x,dmpaB?brst.y,dmpaB?brst.z,dmpaB?.abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','[nT]'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

if 0 % dslE?
  hca = irf_panel('E');
  c_eval('irf_plot(hca,{dslE?.x,dslE?.y},''comp'');',ic)
  hca.YLabel.String = {'E_{DSL}','[mV/m]'};
  irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);
end

hca = irf_panel('brst E');
c_eval('irf_plot(hca,{dslE?brst.x,dslE?brst.y},''comp'');',ic)
hca.YLabel.String = {'E_{DSL}','[mV/m]'};
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst);',ic);
hca.YLabel.String = {'v_e','[km/s]'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca = irf_panel('brst scPot');
c_eval('irf_plot(hca,(-1)*P?brst);',ic);
hca.YLabel.String = '-scPot [V]';

hca = irf_panel('brst n');
c_eval('irf_plot(hca,{ne?_lowres,ni?_lowres},''comp'');',ic);
hca.YLabel.String = {'n','[cm^{-3}]'};
hca.YScale = 'lin';
irf_legend(hca,{'n_e','n_i'},[0.95 0.95]);

if 0 % T?.abs
  hca = irf_panel('brst T');
  c_eval('irf_plot(hca,{Te?_lowres.abs/3,Ti?brst.abs/3},''comp'');',ic);
  hca.YLabel.String = {'T','[eV]'};
  hca.YScale = 'lin';
  irf_legend(hca,{'T_e','T_i'},[0.95 0.95]);
end
if 0
  hca = irf_panel('brst Te');
  c_eval('irf_plot(hca,facTe?);',ic);
  hca.YLabel.String = {'T','[eV]'};
  hca.YScale = 'lin';
  irf_legend(hca,{'T_{\perp,1}','T_{\perp,2}','T_{||}'},[0.95 0.95]);
end

hca = irf_panel('brst Te');
c_eval('irf_plot(hca,Te?_lowres);',ic);
hca.YLabel.String = {'T','[eV]'};
hca.YScale = 'lin';
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.95 0.95]);

if 0
  hca = irf_panel('brst Pe');
  c_eval('irf_plot(hca,facPe?*1e4);',ic);
  hca.YLabel.String = {'P_e','[nPa]'};
  hca.YScale = 'lin';
  irf_legend(hca,{'P_{\perp,1}','P_{\perp,2}','P_{||}'},[0.95 0.95]);
end
if 1
  hca = irf_panel('brst Pe');
  c_eval('irf_plot(hca,Pe?_lowres*1e4);',ic);
  hca.YLabel.String = {'P_e','[nPa]'};
  hca.YScale = 'lin';
  irf_legend(hca,{'P_{x}','P_{y}','P_{z}'},[0.95 0.95]);
end

if 0 % Ti?brst
  hca = irf_panel('brst Ti');
  c_eval('irf_plot(hca,Ti?brst);',ic);
  hca.YLabel.String = 'T_i [eV]';
  hca.YScale = 'lin';
  irf_legend(hca,{'T_x','T_y','T_z'},[0.95 0.95]);
end
if 0 % vi?brst
  hca = irf_panel('brst vi');
  c_eval('irf_plot(hca,vi?brst);',ic);
  hca.YLabel.String = 'v_i [km/s]';
  irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end

hca = irf_panel('Ane');
c_eval('irf_plot(hca,Ane?);',ic);
ylabel(hca,'T_{e,||}/T_{e,\perp}','interpreter','tex');

hca = irf_panel('Agy');
c_eval('irf_plot(hca,Agy?);',ic);
ylabel(hca,'log_{10}(P_{\perp 1}/P_{\perp 2})','interpreter','tex');

if 0
    hca = irf_panel('vExB1');
    irf_plot(hca,VExBav1);
    hca.YLabel.String = 'vi [km/s]';
    irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end

irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h,'y')

%hca = irf_panel('Ane');
%hca.YLim = 100*[-1 1];
%irf_pl_mark(h,tint(1).epochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);
%%
if 0
for ii = 1: numel(h);
    h(ii).Position(3) = h(ii).Position(3)*0.88;
    grid(h(ii),'off')
end
end
%delete(h(end))
