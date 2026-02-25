tint_es = irf.tint('2015-10-16T10:33:29.00Z/2015-10-16T10:33:30.600Z');
c_eval('E?par=dslE?hmfe.dot(dmpaB?brst.resample(dslE?hmfe))/dmpaB?brst.resample(dslE?hmfe).abs;')
c_eval('E?par = irf_filt(E?par.tlim(tint_es),1,0,450,5);')
c_eval('Ve?par=dslE?hmfe.dot(dmpaB?brst.resample(dslE?hmfe))/dmpaB?brst.resample(dslE?hmfe).abs;')
c_eval('Ve?par = irf_filt(E?par.tlim(tint_es),1,0,450,5);')




tint = irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z');
[out,l,v] = irf_minvar(dmpaB1brst.tlim(tint));

L = v(1,:);
M = v(2,:);
N = v(3,:);

c_eval('plve? = irf.ts_vec_xyz(Ve_psd?.time,[Ve_psd?.dot(L).data Ve_psd?.dot(M).data Ve_psd?.dot(N).data]);')

h = irf_plot(2);
if 0
hca=irf_panel('v abs');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{plve1.x,plve2.x,plve3.x,plve4.x},'comp')
hca.YLabel.String = 'v_L (km/s)';
set(hca,'ColorOrder',mms_colors('1234'))
end

hca=irf_panel('ne');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{ne_psd1,ne_psd2,ne_psd3,ne_psd4},'comp')
hca.YLabel.String = 'v_L (km/s)';
set(hca,'ColorOrder',mms_colors('1234'))

hca=irf_panel('E par');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{E1par,E2par,E3par,E4par},'comp')
hca.YLabel.String = 'E_{||} (mV/m)';
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.02 0.95]);

irf_zoom(h,'x',tint_es)
irf_zoom(h,'y')

irf_plot_axis_align
%add_length_on_top(h(1),20,0.5)
%h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'d','e'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0],'fontsize',14)
  %h(ii).Position(3) = h(ii).Position(3)*0.95;
  h(ii).FontSize = 14;
end
hca