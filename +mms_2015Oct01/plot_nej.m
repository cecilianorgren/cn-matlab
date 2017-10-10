
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')
irf_zoom(h,'x',tintOv)
irf_zoom(h,'y')
irf_plot_axis_align
%%

h = irf_plot(8);

hca = irf_panel('n');
irf_plot(hca,{ne1brst,ni1brst},'comp');
hca.YLabel.String = {'n','[cc]'};
irf_legend(hca,{'n_e','n_i'},[0.95 0.95]);

hca = irf_panel('vx');
irf_plot(hca,{ve1brst.x,vi1brst.x},'comp');
hca.YLabel.String = {'v_{y}','[km/s]'};
irf_legend(hca,{'v_e','v_i'},[0.95 0.95]);

hca = irf_panel('vy');
irf_plot(hca,{ve1brst.y,vi1brst.y},'comp');
hca.YLabel.String = {'v_{y}','[km/s]'};
irf_legend(hca,{'v_e','v_i'},[0.95 0.95]);

hca = irf_panel('vz');
irf_plot(hca,{ve1brst.z,vi1brst.z},'comp');
hca.YLabel.String = {'v_{z}','[km/s]'};
irf_legend(hca,{'v_e','v_i'},[0.95 0.95]);

hca = irf_panel('jx');
irf_plot(hca,{je1.x,ji1.x},'comp');
hca.YLabel.String = {'j_{y}','[km/s]'};
irf_legend(hca,{'j_e','j_i'},[0.95 0.95]);

hca = irf_panel('jy');
irf_plot(hca,{je1.y,ji1.y},'comp');
hca.YLabel.String = {'j_{y}','[km/s]'};
irf_legend(hca,{'j_e','j_i'},[0.95 0.95]);

hca = irf_panel('jz');
irf_plot(hca,{je1.z,ji1.z},'comp');
hca.YLabel.String = {'j_{z}','[km/s]'};
irf_legend(hca,{'j_e','j_i'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

if 0
    hca = irf_panel('Jfac');
    irf_plot(hca,jfac);
    %hca.YLabel.String = 'J [nAm^{-2}]';
    ylabel(hca,{'J_{FAC}','[nA m^{-2}]'},'Interpreter','tex');
    irf_legend(hca,{'J_{\perp 1}','J_{\perp 2 (close to y)}','J_{||}'},[0.95 0.95]);
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')

