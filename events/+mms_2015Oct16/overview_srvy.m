tint = irf.tint('2015-10-16T05:30:00.00Z/2015-10-16T16:00:00.00Z');

h = irf_plot(8);

hca = irf_panel('B1');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{dmpaB1.x,dmpaB1.y,dmpaB1.z},'comp');
hca.YLabel.String = {'B','(nT)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'B_x','B_y','B_z',},[0.95 0.95]);

hca = irf_panel('Babs');
irf_plot(hca,dmpaB1.abs);
hca.YLabel.String = {'|B|','(nT)'};
%irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('E1');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{dslE1.x,dslE1.y},'comp');
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);
irf_zoom(hca,'y')

hca = irf_panel('eDEF');
irf_spectrogram(hca,eEnSp1)
hca.YScale = 'log';
hca.YTick = [1e1 1e2 1e3 1e4];
hca.YLim = eEnSp1.f([1 end]);
hca.YLabel.String = {'E_e','(eV)'};

hca = irf_panel('iDEF');
irf_spectrogram(hca,iEnSp1)
%hold(hca,'on');
%irf_plot(hca,vi1.abs);
%irf_plot(hca,VExBav,'w');
hca.YLabel.String = {'E_i','(eV)'};
hca.YScale = 'log';
hca.YTick = [1e1 1e2 1e3 1e4];
hca.YLim = eEnSp1.f([1 end]);
%hold(hca,'off');

hca = irf_panel('n1');
irf_plot(hca,(ni1+ne1)*0.5);
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
hca.YLabel.Interpreter = 'tex';
%irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

hca = irf_panel('vi1');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{vi1.x,vi1.y,vi1.z},'comp');
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'v_i','(km/s)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca = irf_panel('j');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{j.x,j.y,j.z},'comp');
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);


if 0
    hca = irf_panel('vExB1');
    irf_plot(hca,VExBav1);
    hca.YLabel.String = 'vi [km/s]';
    irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end

irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h([1:3 6:8]),'y')
%irf_pl_mark(h,tint(1).epochUnix')
h(1).Title.String = 'MMS 1';
irf_plot_axis_align

for ii = 1: numel(h);
    %h(ii).Position(3) = h(ii).Position(3)*0.88;
    grid(h(ii),'off')
end
%delete(h(end))
%add_position(h(7),Rxyz1)
%% Roy Torberts requested intervals
tint1 = irf.tint('2015-10-16T13:00:00',8*60);
tint2 = irf.tint('2015-10-16T11:25:00',60*7);
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')

%%
tint = irf.tint('2015-10-16T05:00:00.00Z/2015-10-16T16:00:00.00Z'); 
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);',1);

tEpoch = dmpaB1brst.time.epochUnix;
tdiff = diff(tEpoch);
tind = find(tdiff>1);
hhmark = irf_pl_mark(h,dmpaB1brst.time.epochUnix,'green');
%c_eval('dmpaB?brst.name = ''B? brst DMPA'';')
%%
hburstmark=irf_pl_mark(h,[dmpaB1brst.time(tind(1:end-1)+1).epochUnix dmpaB1brst.time(tind(2:end)).epochUnix],'blue');

