w_ratio = 0.1:0.05:1;
mu_diff = 0:50:500;
vt_ratio = 1:0.2:5;

nw = numel(w_ratio);
nmu = numel(mu_diff);
nvt = numel(vt_ratio);

w1 = 1;
vt1 = 100;
mu1 = 0;

d_mc   = zeros(nw,nmu,nvt);
d_grid = zeros(nw,nmu,nvt);
tic
for iw = 1:nw
  w2 = w1*w_ratio(iw);
  for ivt = 1:nvt
    vt2 = vt1*vt_ratio(ivt);
    for imu = 1:nmu
      mu2 = mu1 + mu_diff(imu);
      
      g.w = [w1 w2];
      g.mu = [mu1 0 0; mu2 0 0];
      g.Sigma = cat(3,diag([vt1^2 vt1^2 vt1^2]),diag([vt2^2 vt2^2 vt2^2]));

      [d_mc_, d_grid_] = gaussianity_measure(g,10000,-2000:50:2000);
      d_mc(iw,imu,ivt) = d_mc_;
      d_grid(iw,imu,ivt) = d_grid_;
    end
  end
end
toc

%%
h = setup_subplots(1,3);
isub = 1;

w_range = 1:nw;
mu_range = 1:2;1:nmu;
vt_range = 1:nvt;
if 0
  hca = h(isub); isub = isub + 1;
  x = w_ratio(w_range);
  y = mu_diff(mu_range)/vt1;
  c = squeeze(mean(d_grid,3));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{grid}';
end
if 1
  hca = h(isub); isub = isub + 1;
  x = w_ratio(w_range);
  y = mu_diff(mu_range)/vt1;
  c = squeeze(mean(d_mc(w_range,mu_range,vt_range),3));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{MC}';
end
if 0
  hca = h(isub); isub = isub + 1;
  x = w_ratio(w_range);
  y = vt_ratio(vt_range);
  c = squeeze(mean(d_grid(w_range,mu_range,vt_range),2));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{grid}';
end
if 1
  hca = h(isub); isub = isub + 1;
  x = w_ratio(w_range);
  y = vt_ratio(vt_range);
  c = squeeze(mean(d_mc(w_range,mu_range,vt_range),2));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{MC}';
end
if 0
  hca = h(isub); isub = isub + 1;
  x = mu_diff(mu_range)/vt1;
  y = vt_ratio(vt_range);
  c = squeeze(mean(d_grid(w_range,mu_range,vt_range),1));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{grid}';
end
if 1
  hca = h(isub); isub = isub + 1;
  x = mu_diff(mu_range)/vt1;
  y = vt_ratio(vt_range);
  c = squeeze(mean(d_mc(w_range,mu_range,vt_range),1));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{MC}';
end
c_eval('h(?).Position(2) = h(?).Position(2)+0.05;',1:numel(h))
c_eval('h(?).YDir = "normal";',1:numel(h))
linkprop(h,{'CLim'}); 
h(1).CLim = [0 1];
colormap([flipdim(irf_colormap('Spectral'),1)])
c_eval('h(?).FontSize = 14;',1:numel(h))