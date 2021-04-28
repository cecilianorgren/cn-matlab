%% 
units = irf_units;
mO = units.mp*16; % Oxygen mass

nPolar = 16;
nAzim = 32;
nEnergy = 64;
edgesPolar = linspace(0,180,nPolar+1);
edgesAzim = linspace(0,360,nAzim+1);
edgesEnergy = logspace(log10(10),log10(32000),nEnergy+1); % eV

centerPolar = edgesPolar(1:end-1)+diff(edgesPolar);
centerAzim = edgesAzim(1:end-1)+diff(edgesAzim);
centerEnergy = edgesEnergy(1:end-1)+diff(edgesEnergy);

[EN,AZ,POL] = meshgrid(centerEnergy,centerAzim,centerPolar);
f = EN; % just some random number, whatever your disitrbution happens to be

VABS = sqrt(2*EN*units.eV/mO); % m/s, units.eV = e, tranform from eV to SI units

VX = -VABS.*sin(POL).*cos(AZ); % '-' because the FPI data shows which direction the particles were coming from
VY = -VABS.*sin(POL).*sin(AZ);
VZ = -VABS.*cos(POL);

v_cart = linspace(-fix(max(abs(VX(:)))),fix(max(abs(VX(:)))),20);
[VX_cart,VY_cart,VZ_cart] = meshgrid(v_cart,v_cart,v_cart);
VABS_cart = sqrt(VX_cart.^2 + VY_cart.^2 + VZ_cart.^2);

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

nrows = 1;
ncols = 3;
npanels = nrows*ncols;
isub = 1;

for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

hca = h(isub); isub = isub + 1;
scatter3(hca,EN(1:step:end),AZ(1:step:end),POL(1:step:end),40,VABS(1:step:end),'.');
hca.XLabel.String = 'Energy';
hca.YLabel.String = 'Azimuthal angle';
hca.ZLabel.String = 'Polar angle';
hcb = colorbar('peer',hca);
hcb.YLabel.String = '|v|';

hca = h(isub); isub = isub + 1;
scatter3(hca,VX(1:step:end),VY(1:step:end),VZ(1:step:end),40,VABS(1:step:end),'.'); 
axis(hca,'equal')
hca.XLabel.String = 'v_x';
hca.YLabel.String = 'v_y';
hca.ZLabel.String = 'v_z';
hcb = colorbar('peer',hca);
hcb.YLabel.String = '|v|';

hca = h(isub); isub = isub + 1;
scatter3(hca,VX_cart(1:step:end),VY_cart(1:step:end),VZ_cart(1:step:end),40,VABS_cart(1:step:end),'.'); 
axis(hca,'equal')
hca.XLabel.String = 'v_x';
hca.YLabel.String = 'v_y';
hca.ZLabel.String = 'v_z';
hcb = colorbar('peer',hca);
hcb.YLabel.String = '|v|';

for ip = 1:npanels
  irf_legend(h(ip),labels{ip},[0 1],'fontsize',14)
end

%% Monto carlo vs weighted mean
nv = 20;
vedges = linspace(0,10,nv+1);
vedges = logspace(log10(0.1),log10(10),nv+1);
dv = diff(vedges);
vcenter = vedges(1:end-1)+diff(vedges)/2;
f = 0.1*vcenter+2*rand(1,nv);
f(end) = 0;
n = sum(f.*dv);
%%
% New grid
nv_new = 50;
vedges_new = linspace(0,10,nv_new+1);
dv_new = diff(vedges_new);
vcenter_new = vedges_new(1:end-1)+diff(vedges_new)/2;


% Interpolation
method = 'linear';
f_new_interp = interp1(vcenter,f,vcenter_new,method);
n_new_interp = sum(f_new_interp.*dv_new,'omitnan');

% Monte Carlo particles
n_mc = 3000; % divide 
n_center_new = zeros(1,nv_new);
for iv = 1:nv % loop through bins
  vlow = vedges(iv);
  vhigh = vedges(iv+1);
  v_tmp = vlow + (vhigh-vlow)*rand(n_mc,1);
  dn_tmp = f(iv)*(vhigh-vlow)/n_mc;
  N = histcounts(v_tmp,vedges_new);
  n_center_new = n_center_new + N*dn_tmp;
end
f_new_mc = n_center_new./diff(vedges_new);
n_new_mc = sum(f_new_mc.*dv_new);

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
colors = pic_colors('matlab');

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
isub = 1;

for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
  %h(ip).FontSize = 14;
end

if 0 % original, interpolation
  hca = h(isub); isub = isub + 1;
  f_plot = [f;f];
  v_plot = [vedges(1:end-1);vedges(2:end)];
  hlines = plot(hca,v_plot(:),f_plot(:),vcenter(:),f(:),'*');
  hlines(2).Color = hlines(1).Color;
  %bar(hca,vcenter,f,1,'facealpha',0.5)
  hold(hca,'on')
  hnew = plot(hca,vcenter_new,f_new_interp,'-*');
  %bar(hca,vcenter_new,f_new_interp,'facealpha',0.5)
  hold(hca,'off')
  legend([hlines(1) hnew],{'f_{orig}','f_{interp}'})
end
if 0 % original, monte carlo
  hca = h(isub); isub = isub + 1;
  f_plot = [f;f];
  v_plot = [vedges(1:end-1);vedges(2:end)];
  hlines = plot(hca,v_plot(:),f_plot(:),vcenter(:),f(:),'*');
  hlines(2).Color = hlines(1).Color;
  %bar(hca,vcenter,f,1,'facealpha',0.5)
  hold(hca,'on')
  hnew = plot(hca,vcenter_new,f_new_mc,'-*');
  %bar(hca,vcenter_new,f_new_interp,'facealpha',0.5)
  hold(hca,'off')
  legend([hlines(1) hnew],{'f_{orig}',sprintf('f_{MC}, n_{MC}=%g',n_mc)})
end
if 1 % original, interpolation, monte carlo
  hca = h(isub); isub = isub + 1;
  f_plot = [f;f];
  v_plot = [vedges(1:end-1);vedges(2:end)];
  hlines = plot(hca,v_plot(:),f_plot(:),vcenter(:),f(:),'s');
  hlines(2).Color = hlines(1).Color;
  %bar(hca,vcenter,f,1,'facealpha',0.5)
  hold(hca,'on')
  hinterp = plot(hca,vcenter_new,f_new_interp,'-*');
  hmc = plot(hca,vcenter_new,f_new_mc,'-');
  hmc.Marker = '.';
  %bar(hca,vcenter_new,f_new_interp,'facealpha',0.5)
  hold(hca,'off')
  legend([hlines(1),hinterp,hmc],{sprintf('f_{orig} (%s), n = %.2f',method,n),...
    sprintf('f_{interp}, n = %.2f',n_new_interp),...
    sprintf('f_{MC} (n_{MC}=%g), n = %.2f',n_mc,n_new_mc)},'location','best')
  hca.XLabel.String = 'v';
  hca.YLabel.String = 'f^{1D}';  
end

for ip = 1:npanels  
  h(ip).FontSize = 14;
end