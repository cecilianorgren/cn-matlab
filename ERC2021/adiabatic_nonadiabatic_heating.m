% First adiabatic invariant, 
% mu = mvperp^2/2B, or mu = pperp^2/2mB, where pperp = gamma*m*vperp
% pperp1^2/2mB1 = pperp2^2/2mB2
% sqrt(pperp2_1^2+pperp2_2^2)/sqrt(pperp1_1^2+pperp1_2^2) = sqrt(B2/B1)
% sqrt(2*pperp2_1^2)/sqrt(2*pperp1_1^2) = sqrt(B2/B1)
% pperp2_1/pperp1_1 = sqrt(B2/B1)
% pperp2/pperp1 = sqrt(B2/B1)
% ptot = sqrt(2*pperp^2 + ppar^2)= const;
% ppar = sqrt(ptot^2-2*pperp^2)
% ppar2 = sqrt(ptot^2-2*pperp2^2)
% ppar2/ppar1 = sqrt(ptot^2-2*pperp2^2)/sqrt(ptot^2-2*pperp1^2)

% Second adiabatic invariant, 
% J = int_a^b(ppar)ds, a, b, ends of magnetic bottle

% Beginning point is isotropic Maxwellian distribution
fmax = @(vx,vy,vz,n,vdx,vdy,vdz,vtx,vty,vtz) ...
  n.*(1/pi)^(3/2)/(vtx.*vty.*vtz).*exp(-(vx-vdx).^2./vtx.^2 -(vy-vdy).^2./vty.^2 -(vz-vdz).^2./vtz.^2);

n = 1e-6; % m-3
vt = 1000e3; % m/s
vtpar = vt;
vtperp = vt;
vtx = vtperp;
vty = vtperp;
vtz = vtpar;
vdx = 0; % m/s
vdz = 0; % m/s
vdy = 0; % m/s
nv = 50;
v = 3*vt*linspace(-1,1,nv);

[VX1,VY1,VZ1] = ndgrid(v,v,v);
VTOT = sqrt(VX1.^2 + VY1.^2 + VZ1.^2);
B1 = 10e-9; % T
L = 1000e3; % m/s
F1 = fmax(VX1,VY1,VZ1,n,vdx,vdy,vdz,vtx,vty,vtz);

B2 = 40e-9; % T
VX2 = VX1/sqrt(B2/B1);
VY2 = VY1/sqrt(B2/B1);
VZ2 = VZ1.*sqrt(VTOT.^2 - VX2.^2 - VY2.^2)./sqrt(VTOT.^2 - VX1.^2 - VY1.^2); % ppar2 = sqrt(ptot^2-2*pperp2^2)
F2 = fmax(VX2,VY2,VZ2,n,vdx,vdy,vdz,vtx,vty,vtz);
vx2 = squeeze(VX2(:,1,1));
vy2 = squeeze(VY2(1,:,1));
vz2 = squeeze(VZ2(1,1,:));

nrows = 2;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols); isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,v,v,squeeze(mean(F1,3))');
shading(hca,'flat');
hca.XLabel.String = 'v_{\perp,1}';
hca.YLabel.String = 'v_{\perp,2}';

hca = h(isub); isub = isub + 1;
pcolor(hca,v,v,squeeze(mean(F2,3))');
shading(hca,'flat');
hca.XLabel.String = 'v_{\perp,1}';
hca.YLabel.String = 'v_{\perp,2}';

hca = h(isub); isub = isub + 1;
pcolor(hca,v,v,squeeze(mean(F1,2))');
shading(hca,'flat');
hca.XLabel.String = 'v_{\perp,1}';
hca.YLabel.String = 'v_{||}';

hca = h(isub); isub = isub + 1;
pcolor(hca,v,v,squeeze(mean(F2,2))');
shading(hca,'flat');
hca.XLabel.String = 'v_{\perp,1}';
hca.YLabel.String = 'v_{||}';

hlinks = linkprop(h,{'XLim','YLim','CLim'});

for ipanel = 1:npanels
  hcb = colorbar('peer',h(ipanel));
  axis(h(ipanel),'square')
  axis(h(ipanel),'equal')
  h(ipanel).XLabel.Interpreter = 'tex';
  h(ipanel).YLabel.Interpreter = 'tex';
end


%% Fermi acceleration
wbefore0 = linspace(0,100,100);
v = -1;
wbefore = wbefore0;
N = 3; % number of reflections/energy gains
wmodel = wbefore0*exp(v^2*N);
for iN = 1:N
  vbefore = sqrt(wbefore);
  vafter = -vbefore + 2*v;
  %Wafter = vbefore.^2 - 2*v*vbefore + 4*v.^2; % (-vbefore + 2*v)^2
  % dW = Wafter - Wbefore = vbefore^2 - 2*v*vbefore + 2*v - vbefore.^2 = - 2*v*vbefore + 4*v^2;

  wafter = vafter.^2;
  %dW = -2*v*vbefore + 4*v.^2;
  dw = wafter - wbefore;
  wafter = wbefore + dw;
  
  wbefore = wafter; % for next iteration
end

plot(wbefore0,wafter,wbefore0,wmodel)

