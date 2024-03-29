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
colormap(irf_colormap('waterfall'))

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

%% Map of constant mu(theta1,alpha)
wE = 2*pi/(24*60*60);
RE = 6371200; % m
me = 9.1094e-31;
mp = 1.6726e-27;
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
M = 8.22e22; % A m^2

W_eV = 1e4;
m = me;
v = sqrt(2*W_eV*e/m);



alpha2 = @(B1,B2,alpha1) asind(sqrt(B2./B1).*sind(alpha1));

r = @(r0,theta) r0*sind(theta).^2;
B = @(r,theta) (mu0*M/4/pi./r.^3)*sqrt(3*cosd(theta).^2+1);
B = @(r0,theta) (mu0*M/4/pi./r(r0,theta).^3).*sqrt(3*cosd(theta).^2+1);
alpha = @(r0,theta,alpha0) asind(sqrt(B(r0,theta)./B(r0,90)).*sind(alpha0));

mu = @(r0,theta,alpha) v^2*m*sind(alpha).^2/2./B(r0,theta);

req = 5*RE;
theta = linspace(01,90,300);
alpha0 = linspace(0,90,291);
[AL0,TH] = ndgrid(alpha0,theta);
ALPHA = alpha(req,TH,AL0);
MU = mu(req,TH,ALPHA);

fontsize = 14;
linewidth = 1;
loss_altitude = RE+100e3; % collsision altitude
i_loss_theta = find(r(req,theta)>loss_altitude,1,'first');
loss_theta = theta(i_loss_theta);


nrows = 1; ncols = 2; ipanel = 0;
h = gobjects(nrows,ncols);
for irow = 1:nrows; for icols = 1:ncols; ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end; end 
isub = 1;
if 0 % field line
  hca = h(isub); isub = isub + 1;
  x = r(req,theta).*sind(theta);
  y = r(req,theta).*cosd(theta);
  plot(hca,x/RE,y/RE)
  hold(hca,'on')
  thE = 0:90;
  plot(hca,sind(thE),cosd(thE),'k')
  hold(hca,'off')
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'y (R_E)';
  axis(hca,'equal')
end
if 0 % alpha(r0,theta)
  hca = h(isub); isub = isub + 1; 
  plot(hca,theta,alpha(req,theta,60))  
end
if 0 % B(r0,theta)
  hca = h(isub); isub = isub + 1;
  plot(hca,theta,B(req,theta))  
end
if 1 % r(r0,theta)
  hca = h(isub); isub = isub + 1;
  hl1 = plot(hca,theta,r(req,theta)*1e-3);
  hold(hca,'on')
  hl2 = plot(hca,theta,theta*0+RE*1e-3+100,'k');
  hold(hca,'off')    
  hca.XLabel.String = '\theta (colatitude)';
  hca.XLabel.String = 'Colatitude (\theta, deg)';
  hca.YLabel.String = 'Radial distance (km)';
  hca.XLim = theta([1 end]);
  irf_legend(hca,{sprintf('L = %g', req/RE)},[0.02 0.98],'fontsize',fontsize,'color','k') 
  hold(hca,'on')
  plot(hca,loss_theta*[1 1],hca.YLim,'k--')
  hold(hca,'off')
  legend(hca,{'field line','loss cone altitude','loss cone colatitude'},'box','off','location','best')
end
if 1 % alpha(r0,theta), mu as lines
  hca = h(isub); isub = isub + 1;
  
  plMU = real(MU);
  plMU = MU;
  plMU(abs(imag(ALPHA))>0) = NaN;
  mulev = logspace(-12,-8,10);
  mulev = 20;
  
  plALPHA = real(ALPHA);
  plALPHA = ALPHA;
  plALPHA(abs(imag(plALPHA))>0) = NaN;
  
  pcolor(hca,TH,AL0,plALPHA) 
  shading(hca,'flat')
 
  hold(hca,'on')  
  contour(hca,TH,AL0,plMU,mulev,'k') 
  plot(hca,loss_theta*[1 1],[0 90],'k--')
  hold(hca,'off')
  
  shading(hca,'flat')
  hca.XLabel.String = '\theta (colatitude)';
  hca.XLabel.String = 'Colatitude (\theta, deg)';
  hca.YLabel.String = 'Equatorial pitch angle (\theta_{eq}, deg)';
  hcb = colorbar(hca);
  hcb.Label.String = '\alpha (pitch angle)';
  hcb.Label.String = 'Pitch angle (\alpha, deg)';
  hcb.FontSize = fontsize;
  hca.CLim = [0 90];
  colormap(hca,irf_colormap('waterfall'))
  irf_legend(hca,{sprintf('L = %g', req/RE)},[0.02 0.98],'fontsize',fontsize,'color','k')
  hca.XLim = theta([1 end]);
end
if 0 % mu(r0,theta)
  hca = h(isub); isub = isub + 1;
  plMU = real(MU);
  plMU = MU;
  plMU(abs(imag(ALPHA))>0) = NaN;
  
  mu_scaling = (m/2/e)*1e-6/1e-5; % eV/T -> MeV/G, 
  
  contourf(hca,TH,AL0,plMU/mu_scaling,20) 
  shading(hca,'flat')
  hca.XLabel.String = '\theta (colatitude)';
  hca.YLabel.String = '\alpha_{eq} (equatorial pitch angle)';
  hcb = colorbar(hca);
  hcb.Label.String = '\mu (arb.)';
  hcb.FontSize = fontsize;
  %hca.CLim = [0 90];
  colormap(hca,irf_colormap('waterfall'))
  irf_legend(hca,{sprintf('L = %g', r0/RE)},[0.02 0.98],'fontsize',fontsize,'color','k')
  hca.XLim = theta([1 end]);
  hold(hca,'on')
  plot(hca,loss_theta*[1 1],[0 90],'k--')
  hold(hca,'off')
  hcb.YTick = [];
  %hcb.Location = 'north';
  %hcb.Position = [0.4 0.27 0.3 0.02];
end
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).Position(2) = 0.20;',1:numel(h))
c_eval('h(?).Position(4) = 0.7;',1:numel(h))
%irf_plot_axis_align(h)