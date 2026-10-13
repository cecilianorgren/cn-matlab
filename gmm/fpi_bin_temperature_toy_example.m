% 
% 
% 
%% Compare energy resolution to velocity resolution
% Generate energy grid
eE = logspace(1,log10(32e3),33);

E1 = eE(1:end-1); % lower edges of bins
E2 = eE(2:end); % upper edges of bin
dE = diff(eE); % bin widths
E0 = eE(1:end-1) + dE; % center energies of bins

% Generate velocity grid from energy grid
eV = sqrt(2*eE); % edges
V0 = sqrt(2*E0); % centers, corresponding to the energy centers
V1 = sqrt(2*E1); % lower edges
V2 = sqrt(2*E2); % upper edges
dV = V2-V1; % bin width


Eres = dE./E0; % energy resolution
Vres = dV./V0; % speed resolution (recall that v0 is not center of speed grid)

% plot(dVcV./dEcE)

%% Check temperature of bin
Tfactor = 2; % kBT = mvt^2/Tfactor, definition of temperature/thermal speed

% Construct grid
E0 = 1000; % eV
Eres = 0.24;
dE = E0*Eres;
E1 = E0-0.5*dE;
E2 = E0+0.5*dE;

% v0E = sqrt(E0*2*units.eV/units.mp)*1e-3;
% v1 = sqrt(E1*2*units.eV/units.mp)*1e-3;
% v2 = sqrt(E2*2*units.eV/units.mp)*1e-3;
v0E = sqrt(E0*2);
v1 = sqrt(E1*2);
v2 = sqrt(E2*2);


dv = v2-v1;
v0 = v1 + 0.5*dv;
vres = dv/v0;
disp("-------------------------")
disp("E0 = " + E0 + ", Eres = " + Eres + ", dE = Eres*E0 = " + dE)
disp("dv_bin = " + dv + ", vres = " + vres + ", vres/Eres = " + vres/Eres)

% Generate particles
rng(1)
N = 1000000;
p = 1;
r = rand(N,1);
v = (v1.^p + r*(v2.^p-v1.^p)).^(1/p);
vvar = var(v); % sum((v-vmean)^2), T is propto var
vt = sqrt(vvar)/Tfactor;
vmean = mean(v);
disp("var(v) = " + vvar + ", vt = sqrt(var(v))/" + Tfactor + " = " + vt + ", dv_bin = " + dv + ", (dv_bin)^2/12 = " + dv^2/12 + " % Note: visually, 'what one compares' is vt and dv_bin/2")

% T ~ 0.5*vvar ((v-vmean)^2/2)
% Tbin = vvar*1e6*units.mp/units.eV/Tfactor;
Tbin = vvar/Tfactor;
disp("Tbin = var(v)/2 = " + Tbin + ", Ebin = " + dE + ", Tbin/dEbin = " + Tbin/dE)

% Compare to the variance a gaussian fit gives
g = fitgmdist(v,1);
g.Sigma; % ~vt^2
%g.Sigma/2/dE
disp("Sigma (=vt^2) from fitgmdist = " + g.Sigma)

% Plot
hca = subplot(1,1,1);
vgrid = linspace(v1-0.5*dv,v2+0.5*dv,100);
histogram(hca,v,vgrid,'Normalization','probability')

% Compare to normal distribution with the same variance
gaussianPdf = normpdf(vgrid, vmean, sqrt(vvar));
% Normalize the Gaussian curve to the histogram bin width
binWidth = vgrid(2) - vgrid(1);
gaussianPdf = gaussianPdf * binWidth;

vn = vmean + randn(N,1)*sqrt(vvar);
hold(hca,'on')
plot(hca,vgrid,  gaussianPdf, 'r-', 'LineWidth', 1.5)
histogram(hca,vn,vgrid,'Normalization','probability')
patch(hca,vmean + vt*[-1 1 1 -1], hca.YLim(2)*[0 0 1 1],'g','facealpha',0.2)
hca.Children = circshift(hca.Children,-1);

legend(hca,'-vt < v - vmean < vt','Box function samples', 'Equal-variance Gaussian', 'Equal-variance Gaussian samples','location','northoutside')
%hca.Title.String = "E_0 = " + E0 + "eV, dE = " + dE + ", dE/E0 = " + Eres + ", k_BT = mv^2/" + Tfactor;
hca.Title.String = "E_0 = " + E0 + ", dE = " + dE + ", dE/E0 = " + Eres + ", T = v_t^2/" + Tfactor;
hold(hca,'off')
