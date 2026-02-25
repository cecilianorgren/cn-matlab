% Alfven phase velocity and groupd velocity

VA = 4000e3;
CS = 1000e3;
CS = VA/10;

if 1
  kx = linspace(0,1,100); dkx = kx(2)-kx(1);
  ky = linspace(0,1,100); dky = ky(2)-ky(1);
  [KX,KY] = meshgrid(kx,ky);
  K = sqrt(KX.^2 + KY.^2);
  THETA = atand(KY./KX);
else
  k = linspace(0,1,20); dk = k(2)-k(1);
  theta = linspace(0,360,20);

  [K,THETA] = meshgrid(k,theta);
  KX = K.*cosd(THETA);
  KY = K.*sind(THETA);
end

w1 = @(k,theta) sqrt(k.^2*VA.*cosd(theta).^2);
w2 = @(k,theta) sqrt(0.5*k.^2.*(VA^2+CS^2+sqrt((VA^2+CS^2)^2-4*VA^2*CS^2.*cosd(theta).^2)));
w3 = @(k,theta) sqrt(0.5*k.^2.*(VA^2+CS^2-sqrt((VA^2+CS^2)^2-4*VA^2*CS^2.*cosd(theta).^2)));

vp1 = @(k,theta) sqrt(VA.*cosd(theta).^2);
vp2 = @(k,theta) sqrt(0.5*(VA^2+CS^2+sqrt((VA^2+CS^2)^2-4*VA^2*CS^2.*cosd(theta).^2)));
vp3 = @(k,theta) sqrt(0.5*(VA^2+CS^2-sqrt((VA^2+CS^2)^2-4*VA^2*CS^2.*cosd(theta).^2)));

W1 = w1(K,THETA);
W2 = w2(K,THETA);
W3 = w3(K,THETA);

VP1 = vp1(K,THETA);
VP2 = vp2(K,THETA);
VP3 = vp3(K,THETA);

[vgx1,vgy1] = gradient(W1,KX,KY);
[vgx2,vgy2] = gradient(W2,KX,KY);
[vgx3,vgy3] = gradient(W3,KX,KY);

[vgx1,vgy1] = gradient(W1,K,THETA);
[vgx2,vgy2] = gradient(W2,K,THETA);
[vgx3,vgy3] = gradient(W3,K,THETA);

[vgx1,vgy1] = gradient(W1,kx,ky);
[vgx2,vgy2] = gradient(W2,kx,ky);
[vgx3,vgy3] = gradient(W3,kx,ky);
q

% Plot for exercise
nrows = 1;
ncols = 3;
npanels = nrows*ncols;
for ipl = 1:npanels
  h(ipl) = subplot(nrows,ncols,ipl);
end
isub = 1;

qtpx = 4:10:numel(kx); % quivers to plot
qtpy = 4:10:numel(ky);
hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W1); hold(hca,'on')
quiver(hca,KX(qtpx,qtpy),KY(qtpx,qtpy),vgx1(qtpx,qtpy),vgy1(qtpx,qtpy)); hold(hca,'off')

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W2); hold(hca,'on')
quiver(hca,KX(qtpx,qtpy),KY(qtpx,qtpy),vgx2(qtpx,qtpy),vgy2(qtpx,qtpy)); hold(hca,'off')

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W3); hold(hca,'on')
quiver(hca,KX(qtpx,qtpy),KY(qtpx,qtpy),vgx3(qtpx,qtpy),vgy3(qtpx,qtpy)); hold(hca,'off')

for ip = 1:npanels
  axis(h(ip),'equal')
  axis(h(ip),'square')
  h(ip).XLabel.String = 'k_{||}';
  h(ip).YLabel.String = 'k_{\perp}';
  h(ip).XLim = [0 1];
  h(ip).YLim = [0 1];
  h(ip).FontSize = 16;
  h(ip).Title.String = irf_ssub('\omega_?(k)',ip);
  %h(ip).CLim = h(1).CLim;
end
ip = 2;
h(ip).Title.String = {sprintf('v_A/c_s = %.2f',VA/CS),h(ip).Title.String};

%%
nrows = 2;
ncols = 3;
npanels = nrows*ncols;
for ipl = 1:npanels
  h(ipl) = subplot(nrows,ncols,ipl);
end
isub = 1;

if 0 % w(k) surfaces
  hca = h(isub); isub = isub + 1;
  surf(hca,KX,KY,W1)
  hca.XLabel.String = 'k_x';
  hca.YLabel.String = 'k_y';

  hca = h(isub); isub = isub + 1;
  surf(hca,KX,KY,W2)
  hca.XLabel.String = 'k_x';
  hca.YLabel.String = 'k_y';

  hca = h(isub); isub = isub + 1;
  surf(hca,KX,KY,W3)
  hca.XLabel.String = 'k_x';
  hca.YLabel.String = 'k_y';
end

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,VP1); hold(hca,'on')
quiver(hca,KX,KY,vgx1,vgy1); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,VP2); hold(hca,'on')
quiver(hca,KX,KY,vgx2,vgy2); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,VP3); hold(hca,'on')
quiver(hca,KX,KY,vgx3,vgy3); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W1); hold(hca,'on')
quiver(hca,KX,KY,vgx1,vgy1); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W2); hold(hca,'on')
quiver(hca,KX,KY,vgx2,vgy2); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

hca = h(isub); isub = isub + 1;
contour(hca,KX,KY,W3); hold(hca,'on')
quiver(hca,KX,KY,vgx3,vgy3); hold(hca,'off')
hca.XLabel.String = 'k_x';
hca.YLabel.String = 'k_y';

for ip = 4:npanels
  axis(h(ip),'equal')
end
