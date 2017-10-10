% Plot different electric field products, to compare them and see if they
% agree or not. 


% Plots single satellite overview of diffusion region using burst data
% First run mms_Oct16.load_brst_data

tint = irf.tint('2015-10-16T10:33:10.00Z/2015-10-16T10:33:50.00Z');
ic = 1;


if 1 % make electric field wavelet
  tic
  c_eval('wavE? = irf_wavelet(dslE?brst.abs);',ic)
  toc
end

h = irf_plot(9);
tic
hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?.x,dmpaB?.y,dmpaB?.z,dmpaB?.abs},''comp'');',ic)
hca.YLabel.String = 'B_{DMPA} [nT]';
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('E');
c_eval('irf_plot(hca,{dslE?.x,dslE?.y},''comp'');',ic)
hca.YLabel.String = 'E_{DSL} [mV/m]';
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

hca = irf_panel('brst E');
c_eval('irf_plot(hca,{dslE?brst.x,dslE?brst.y},''comp'');',ic)
hca.YLabel.String = 'E_{DSL} [mV/m]';
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

hca = irf_panel('brst vx');
c_eval('irf_plot(hca,{vExB?brst.x,vi?brst.x,ve?brst.x},''comp'');',ic);
hca.YLabel.String = 'v_x [km/s]';
irf_legend(hca,{'v_{ExB}','v_i','v_e'},[0.95 0.95]);
hca = irf_panel('brst vy');
c_eval('irf_plot(hca,{vExB?brst.y,vi?brst.y,ve?brst.y},''comp'');',ic);
hca.YLabel.String = 'v_y [km/s]';
irf_legend(hca,{'v_{ExB}','v_i','v_e'},[0.95 0.95]);
hca = irf_panel('brst vz');
c_eval('irf_plot(hca,{vExB?brst.z,vi?brst.z,ve?brst.z},''comp'');',ic);
hca.YLabel.String = 'v_z [km/s]';
irf_legend(hca,{'v_{ExB}','v_i','v_e'},[0.95 0.95]);

hca = irf_panel('wavelet E');
c_eval('irf_spectrogram(hca,wavE);',ic);
hold(hca,'on');
c_eval('irf_plot(hca,flh?*1e-3);',ic);
c_eval('irf_plot(hca,fce?*1e-3);',ic);
c_eval('irf_plot(hca,fcp?*1e-3);',ic);
c_eval('irf_plot(hca,fpp?*1e-3);',ic);
c_eval('irf_plot(hca,fpe?*1e-3);',ic);
hold(hca,'off');

hca.YLabel.String = 'f (kHz)';
hca.YScale = 'log';
hca.YLim = [wavE.f(end) wavE.f(1)]*1e-3;
hca.YTick = [1e1 1e2 1e3 1e4]*1e-3;

irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h,'y')
%irf_pl_mark(h,tint(1).epochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);
toc
hca = irf_panel('wavelet E');
hca.YLim = [wavE.f(end) wavE.f(1)]*1e-3;
%%

hca = irf_panel('brst vi');
c_eval('irf_plot(hca,vi?brst);',ic);
hca.YLabel.String = 'v_i [km/s]';
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst);',ic);
hca.YLabel.String = 'v_e [km/s]';
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

if 0
    hca = irf_panel('vExB1');
    irf_plot(hca,VExBav1);
    hca.YLabel.String = 'vi [km/s]';
    irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end

irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h,'y')
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
