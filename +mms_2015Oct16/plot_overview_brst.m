% Plots single satellite overview of diffusion region using burst data
% First run mms_Oct16.load_brst_data

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z');
ic = 1;

h = irf_plot(10);

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

hca = irf_panel('brst scPot');
c_eval('irf_plot(hca,(-1)*P?brst);',ic);
hca.YLabel.String = '-scPot [V]';

if 0
hca = irf_panel('eDEF');
irf_spectrogram(hca,eEnSp1)
hca.YScale = 'log';
hca.YTick = [1e1 1e2 1e3 1e4];
hca.YLim = eEnSp1.f([1 end]);
hca.YLabel.String = 'E_e [eV]';
end
if 0
hca = irf_panel('iDEF');
irf_spectrogram(hca,iEnSp1)
%hold(hca,'on');
%irf_plot(hca,vi1.abs);
%irf_plot(hca,VExBav,'w');
hca.YLabel.String = 'E_i [eV]';
hca.YScale = 'log';
hca.YTick = [1e1 1e2 1e3 1e4];
hca.YLim = eEnSp1.f([1 end]);
%hold(hca,'off');
end
hca = irf_panel('brst n');
c_eval('irf_plot(hca,{ne?brst,ni?brst},''comp'');',ic);
hca.YLabel.String = 'n [cm^{-3}]';
hca.YScale = 'lin';
irf_legend(hca,{'n_e','n_i'},[0.95 0.95]);



hca = irf_panel('brst T');
c_eval('irf_plot(hca,{Te?brst.abs/3,Ti?brst.abs/3},''comp'');',ic);
hca.YLabel.String = 'T [eV]';
hca.YScale = 'lin';
irf_legend(hca,{'T_e','T_i'},[0.95 0.95]);

hca = irf_panel('brst Te');
c_eval('irf_plot(hca,Te?brst);',ic);
hca.YLabel.String = 'T_e [eV]';
hca.YScale = 'lin';
irf_legend(hca,{'T_x','T_y','T_z'},[0.95 0.95]);

hca = irf_panel('brst Ti');
c_eval('irf_plot(hca,Ti?brst);',ic);
hca.YLabel.String = 'T_i [eV]';
hca.YScale = 'lin';
irf_legend(hca,{'T_x','T_y','T_z'},[0.95 0.95]);

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