w_ratio = 0.1:0.05:1;
mu_diff = 0:100:500;
vt_ratio = 1:0.2:5;

w_ratio = logspace(-1.5,1,10); w_scale = 'log';
mu_diff = linspace(0,2,11); mu_scale = 'lin';
vt_ratio = logspace(-1,0,12); vt_scale = 'log';


nw = numel(w_ratio);
nmu = numel(mu_diff);
nvt = numel(vt_ratio);

w1 = 1;
mu1 = 0;
vt1 = 1;

d_mc   = zeros(nw,nmu,nvt);
d_grid = zeros(nw,nmu,nvt);
vgrid = linspace(-3*vt1,3*vt1,100);
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

      [d_mc_, d_grid_] = gaussianity_measure(g,10000,vgrid);
      d_mc(iw,imu,ivt) = d_mc_;
      d_grid(iw,imu,ivt) = d_grid_;
    end
  end
end
toc

%%
h = setup_subplots(1,3);
isub = 1;

w_range = 3;%1:nw;
mu_range = 6;%1:nmu;
vt_range = 6;%1:nvt;

if 1 % grid
  hca = h(isub); isub = isub + 1;
  x = w_ratio(:);
  y = mu_diff(:)/vt1;
  c = squeeze(mean(d_grid(:,:,vt_range),3));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'errL1';
  hca.XScale = w_scale;
  hca.YScale = mu_scale;
  hca.Title.String = "v_{t2}/v_{t1} = " + num2str(vt_ratio(vt_range));
end
if 0 % MC
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
if 1 % grid
  hca = h(isub); isub = isub + 1;
  x = w_ratio(:);
  y = vt_ratio(:);
  c = squeeze(mean(d_grid(:,mu_range,:),2));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'errL1';
  hca.XScale = w_scale;
  hca.YScale = vt_scale;
  hca.Title.String = "(\mu_2-\mu_1)/v_{t1} = " + num2str(mu_diff(mu_range));
end
if 0 % MC
  hca = h(isub); isub = isub + 1;
  x = w_ratio(:);
  y = vt_ratio(:);
  c = squeeze(mean(d_mc(w_range,:,:),1));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = 'w_2/w_1';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'G_{MC}';
end
if 1 % grid
  hca = h(isub); isub = isub + 1;
  x = mu_diff(:)/vt1;
  y = vt_ratio(:);
  c = squeeze(mean(d_grid(w_range,:,:),1));
  pcolor(hca,x,y,c')
  shading(hca,'flat')
  hca.XLabel.String = '(\mu_2-\mu_1)/v_{t1}';
  hca.YLabel.String = 'v_{t2}/v_{t1}';
  hcb = colorbar(hca);
  hcb.Label.String = 'errL1';
  hca.XScale = mu_scale;
  hca.YScale = vt_scale;
  hca.Title.String = "w_{2}/w_{1} = " + num2str(w_ratio(w_range));
end
if 0 % mc
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
c_eval('h(?).Position(4) = h(?).Position(4)-0.05;',1:numel(h))
%c_eval('h(?).YDir = "normal";',1:numel(h))
linkprop(h,{'CLim'}); 
h(1).CLim = [0 .3];
colormap([flipdim(irf_colormap('Spectral'),1)])
c_eval('h(?).FontSize = 10;',1:numel(h))

%% contours at given values
h = setup_subplots(3,1);
isub = 1;

lev = 0.2; 
w_range = 1:nw;
mu_range = 1:nmu;
vt_range = 1:nvt;
if 1 % grid
  colors  = irf_colormap('spectral','interp',nvt);
  hca = h(isub); isub = isub + 1;
  for i = 1:nvt    
    x = w_ratio(w_range);
    y = mu_diff(mu_range)/vt1;
    c = squeeze(d_grid(:,:,i));
    contour(hca,x,y,c',lev*[1 1],'color',colors(i,:))
    hca.XLabel.String = 'w_2/w_1';
    hca.YLabel.String = '(\mu_2-\mu_1)/v_{t1}';
    if i == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
end
if 1 % grid
  colors  = irf_colormap('spectral','interp',nmu);
  hca = h(isub); isub = isub + 1;
  for i = 1:nmu
    x = w_ratio(w_range);
    y = vt_ratio(vt_range);
    c = squeeze(d_grid(:,i,:));
    contour(hca,x,y,c',lev*[1 1],'color',colors(i,:))
    hca.XLabel.String = 'w_2/w_1';
    hca.YLabel.String = 'v_{t2}/v_{t1}';
    if i == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
end
if 1 % grid
  colors  = irf_colormap('spectral','interp',nw);
  hca = h(isub); isub = isub + 1;
  for i = 1:nw
    x = mu_diff(mu_range);
    y = vt_ratio(vt_range);
    c = squeeze(d_grid(i,:,:));
    contour(hca,x,y,c',lev*[1 1],'color',colors(i,:))
    hca.YLabel.String = 'v_2/v_1';
    hca.XLabel.String = '(\mu_2-\mu_1)/v_{t1}';
    if i == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
end