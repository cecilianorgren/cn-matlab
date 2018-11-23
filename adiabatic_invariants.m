% adiabatic_invariants

% Isotropic distribution
np = 1000;
angle = 0  + (180-0)*rand(np,1);
nv = 100;
vmax = 1;
v = 0  + (vmax-0)*rand(np,1);
E = v.^2/2;
phi = 0.03;

vperp0 = v.*sind(angle);
vpar0 = v.*cosd(angle);
B0 = 1;

% mu*B = Eperp;
mu = vperp0.^2/2/B0;

vperp = @(B) sqrt(2*mu*B);
vpar = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);
vpar_trapped = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);
vpar_passing = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);


B1 = 0.5*B0;

figure(200)
nrows = 2;
ncols = 2;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

hca = h(isub); isub = isub + 1;
scatter(hca,vpar0,vperp0)
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vpar(B1),vperp(B1))
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vpar0,vpar(B1))
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vperp0,vperp(B1))
axis(hca,'equal')
%hca.YLim = [0 vmax];
%hca.XLim = [-vmax vmax];


%% Egedal 2008, p. 9-10-
% distribution function
ninfty = 1;
vtpar = 1;
vtperp = 1;
B0 = 1;
Bloc = 0.5;
phi = 1;

f = @(vpar,vperp1,vperp2) ninfty/((pi)^(3/2)*vtpar*vtperp^2)*exp(-(vpar.^2/vtpar^2+vperp1.^2/vtperp^2+vperp2.^2/vtperp^2));

%f = @(vpar,vperp1,vperp2) ninfty/((pi)^(3/2)*vtpar*vtperp^2)* exp(-vperp.^2/vtpar^2*(1-B/Bloc)-phi/vtpar^2-(B0/Bloc)*(vperp1.^2/vtperp^2+vperp2.^2/vtperp^2));

nvperp = 110; 
nvpar = 110;
vmax = *vtpar;
vperp = linspace(-vmax,vmax,nvperp);
vperp1 = vperp;
vperp2 = vperp;
vpar = linspace(-vmax,vmax,nvpar);
dvperp = vperp1(2)-vperp1(1);
dvpar = vpar(2)-vpar(1);
d3v = dvpar*dvperp^2; 

[VPAR,VPERP1,VPERP2] = meshgrid(vpar,vperp1,vperp2);
F = f(VPAR,VPERP1,VPERP2);
nint = sum(F(:))*d3v;

hca = subplot(2,2,1);
pcolor(hca,squeeze(VPAR(:,:,1))-0,squeeze(VPERP1(:,:,1)),sum(F,3));
%hca.Title.String = sprintf('ninfty = %g, nint = %g',ninfty,nint);
hca.Title.String = sprintf('original');
hca.XLim = vmax*[-1 1];
hca.YLim = vmax*[-1 1];
shading(hca,'flat')
hca.YLabel.String = 'v_\perp';
hca.XLabel.String = 'v_{||}';
colorbar;

hca = subplot(2,2,2);
pcolor(hca,squeeze(VPAR(:,:,1))-1,squeeze(VPERP1(:,:,1)),sum(F,3));
%hca.Title.String = sprintf('ninfty = %g, nint = %g',ninfty,nint);
hca.Title.String = sprintf('accelerated');
hca.XLim = vmax*[-1 1];
hca.YLim = vmax*[-1 1];
shading(hca,'flat')
hca.YLabel.String = 'v_\perp';
hca.XLabel.String = 'v_{||}';
colorbar;

[VPAR,VPERP1,VPERP2] = meshgrid(vpar,vperp1,vperp2);
F = f(VPAR*0.5,VPERP1*0.5,VPERP2*0.5);
nint = sum(F(:))*d3v;

hca = subplot(2,2,3);
pcolor(hca,squeeze(VPAR(:,:,1))-0,squeeze(VPERP1(:,:,1)),sum(F,3));
%hca.Title.String = sprintf('ninfty = %g, nint = %g',ninfty,nint);
hca.Title.String = sprintf('heated');
hca.XLim = vmax*[-1 1];
hca.YLim = vmax*[-1 1];
shading(hca,'flat')
hca.YLabel.String = 'v_\perp';
hca.XLabel.String = 'v_{||}';
colorbar;

hca = subplot(2,2,4);
pcolor(hca,squeeze(VPAR(:,:,1))-1,squeeze(VPERP1(:,:,1)),sum(F,3));
%hca.Title.String = sprintf('ninfty = %g, nint = %g',ninfty,nint);
hca.Title.String = sprintf('accelerated + heated');
hca.XLim = vmax*[-1 1];
hca.YLim = vmax*[-1 1];
shading(hca,'flat')
hca.YLabel.String = 'v_\perp';
hca.XLabel.String = 'v_{||}';
colorbar;
%axis square
%axis equal