%% 2D mat of theta2
theta2 = @(B1,B2,theta1) asind(sqrt(B2./B1).*sind(theta1));

theta1 = 1:1.5:90;
theta1 = logspace(-1,log10(90),1000);
b1 = 1;
b2 = logspace(-1,log10(3000),2000);
b2b1 = b2./b1;

[B2B1,THETA1] = ndgrid(b2b1,theta1);
THETA2 = real(theta2(b1,b1*B2B1,THETA1));
THETA2(THETA2==90) = NaN;

nrows = 1; ncols = 2; ipanel = 0;
h = gobjects(nrows,ncols);
for irow = 1:nrows; for icols = 1:ncols; ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end; end 
isub = 1;
if 1 % linear scales
  hca = h(isub); isub = isub + 1;
  pcolor(hca,B2B1,THETA1,THETA2)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.Label.String = '\theta_2 (deg.)';
  hca.XLabel.String = 'B_2/B_1';
  hca.YLabel.String = '\theta_1 (deg.)';
  hca.FontSize = 14;
  hca.Color = 0.9*[1 1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XTick = 0:500:10000;
  hca.XTickLabelRotation = 0;
end
if 1 % log B  x-scale
  hca = h(isub); isub = isub + 1;
  pcolor(hca,B2B1,THETA1,THETA2)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.Label.String = '\theta_2 (deg.)';
  hca.XLabel.String = 'B_2/B_1';
  hca.YLabel.String = '\theta_1 (deg.)';
  hca.FontSize = 14;
  hca.Color = 0.9*[1 1 1];
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XTick = 10.^[-1:5];
end
c_eval('h(?).Position([2 4]) = [0.20 0.75];',1:numel(h))
%% Loss cone as a function of B1/B1
theta1 = @(B1,B2,theta2) asind(sqrt(B1./B2).*sind(theta2));

losscone = @(B2B1,theta2) asind(sqrt(1./B2B1).*sind(theta2));

b1 = 1;
b2 = logspace(-1,log10(3000),2000);
b2b1 = b2./b1;

hca = subplot(1,1,1);
plot(hca,b2b1,losscone(b2b1,90))


%% Line plot of theta2
theta2 = @(B1,B2,theta1) asind(sqrt(B2./B1).*sind(theta1));

theta1 = 0:90;
b1 = 1;
b2 = logspace(-1,3,100);
b2b1 = b2./b1;

hca = subplot(1,1,1);
plot(hca,b2b1,theta2(b1,b1*b2b1,theta1))

%% 2D mat of B2(B1,theta1,theta2=90)
theta2 = @(B1,B2,theta1) asind(sqrt(B2./B1).*sind(theta1));

theta1 = 1:1.5:90;
theta1 = logspace(-1,log10(90),1000);
b1 = 1;
b2 = logspace(-1,log10(3000),2000);
b2b1 = b2./b1;

[B2B1,THETA1] = ndgrid(b2b1,theta1);
THETA2 = real(theta2(b1,b1*B2B1,THETA1));
THETA2(THETA2==90) = NaN;

nrows = 1; ncols = 2; ipanel = 0;
h = gobjects(nrows,ncols);
for irow = 1:nrows; for icols = 1:ncols; ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end; end 
isub = 1;
if 1 % linear scales
  hca = h(isub); isub = isub + 1;
  pcolor(hca,B2B1,THETA1,THETA2)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.Label.String = '\theta_2 (deg.)';
  hca.XLabel.String = 'B_2/B_1';
  hca.YLabel.String = '\theta_1 (deg.)';
  hca.FontSize = 14;
  hca.Color = 0.9*[1 1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XTick = 0:500:10000;
  hca.XTickLabelRotation = 0;
end
if 1 % log B  x-scale
  hca = h(isub); isub = isub + 1;
  pcolor(hca,B2B1,THETA1,THETA2)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.Label.String = '\theta_2 (deg.)';
  hca.XLabel.String = 'B_2/B_1';
  hca.YLabel.String = '\theta_1 (deg.)';
  hca.FontSize = 14;
  hca.Color = 0.9*[1 1 1];
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XTick = 10.^[-1:5];
end
c_eval('h(?).Position([2 4]) = [0.20 0.75];',1:numel(h))
