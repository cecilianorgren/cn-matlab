clear dt_ind ind_diff v

tintBL = irf.tint('2015-11-12T07:19:20.70Z/2015-11-12T07:19:21.50');
c_eval('B?_ = mvaB?.resample(mvaE1).tlim(tintBL); B?_mat = B?_.data;',1:4)
timeline = B1_.time;

nB = 100;
Bstep = linspace(-9,8,nB);

dt_ind = nan(nB,4);

for iB = 1:nB
  c_eval('ind? = find(abs(B?_mat(:,1)-Bstep(iB)) == min(abs(B?_mat(:,1)-Bstep(iB)))); dt_ind(iB,?) = ind?;',1:4)  
end

c_eval('ind_diff(:,?) = dt_ind(:,?)-dt_ind(:,1);',1:4)
dt = timeline(2) - timeline(1);
dt_diff = ind_diff*dt;
timeline_diff = timeline(dt_ind(:,1));
ts_dt_diff = irf.ts_scalar(timeline_diff,dt_diff);
  
% get v
v = zeros(nB,3);
for it = 1:nB
  %v(it,:) = irf_4_v([0 mvaRR1],[0 mvaRR2],[0 mvaRR3],[0 mvaRR4],dt_diff(it,:)+10*299792.458);
  v(it,:) = irf_4_v(mvaRR1,mvaRR2,mvaRR3,mvaRR4,dt_diff(it,:)+10*299792.458);
end

ts_v = irf.ts_scalar(timeline_diff,v);

% Plot
npanels = 3;
h = irf_plot(npanels);

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_plot(hca,{mvaB1.x.tlim(tintBL),mvaB2.x.tlim(tintBL),mvaB3.x.tlim(tintBL),mvaB4.x.tlim(tintBL)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % dt
  hca = irf_panel('dt');
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_plot(hca,ts_dt_diff);
  hca.YLabel.String = {'dt','(s)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % v
  hca = irf_panel('v');
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_plot(hca,ts_v);
  hca.YLabel.String = {'v','(s)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end

irf_zoom(h,'x',tintBL)
irf_zoom(h,'y')
irf_plot_axis_align