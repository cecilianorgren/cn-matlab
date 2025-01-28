tint = irf.tint('2015-10-16T10:33:47.00Z/2015-10-16T10:33:50.00Z'); % magnetosphere-magnetosheath-magnetosphere
dtBz = [0 0.03 -0.04 0.06]; vBz = 145*[-0.50 0.4 -0.7];
dtBy = [0 0 -0.13 -0.02];
dtBx = [0 0.04 -0.03 0.07];

dt = dtBz;
v = vBz;

h = irf_plot(4);
hca = irf_panel('will delete'); 


hca = irf_panel('Bx');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.x.tlim(tint),dmpaB2brst.x.tlim(tint),dmpaB3brst.x.tlim(tint),dmpaB4brst.x.tlim(tint)},'comp','dt',dt);
hca.YLabel.String = {'B_{x}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[0.01 0.9]);
irf_legend(hca,{['dt = ' num2str(dt(1))],[' ' num2str(dt(2))],[' ' num2str(dt(3))],[' ' num2str(dt(4))]},[0.01 0.75]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.y.tlim(tint),dmpaB2brst.y.tlim(tint),dmpaB3brst.y.tlim(tint),dmpaB4brst.y.tlim(tint)},'comp','dt',dt);
hca.YLabel.String = {'B_{y}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
%irf_legend(hca,{['dt = ' num2str(dt(1))],[' ' num2str(dt(2))],[' ' num2str(dt(3))],[' ' num2str(dt(4))]},[1.01 0.7]);

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2brst.z.tlim(tint),dmpaB3brst.z.tlim(tint),dmpaB4brst.z.tlim(tint)},'comp','dt',dt);
hca.YLabel.String = {'B_{z}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
%irf_legend(hca,{['dt = ' num2str(dt(1))],[' ' num2str(dt(2))],[' ' num2str(dt(3))],[' ' num2str(dt(4))]},[1.01 0.7]);

irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
for ii = 1:numel(h); 
%  h(ii).Position(3) = h(ii).Position(3)*0.8; 
  h(ii).GridColor = [1 1 1]; 
end

delete(h(1));  

ax = add_length_on_top(h(2),norm(v),0.5);
ax.XLabel.String = 'km';
ax.Title.String = {['velocity = ' num2str(norm(v),'%.0f') '\times[ ' num2str(v(1)/norm(v),'%.2f') ' ' num2str(v(2)/norm(v),'%.2f') ' ' num2str(v(3)/norm(v),'%.2f')  '] km/s']};
