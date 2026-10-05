%posterior = gm{end}.posterior(R);

iK = 2;
[tsN,tsV,tsT,tsFx,tsFy,tsFz,tsN_tot,tsFx_tot,tsFy_tot,tsFz_tot] = gmm_get_moments_and_dists(times,gm(:,iK),Fobs.mid{:},ntot,'sort',true);


%Ek = cellfun(@(x) units.mp*sum(x.mu.^2,2)/2/units.eV,gm);
Ek = cellfun(@(x) units.mp*sum(x.mu.^2*1e6,2)/2/units.eV,gm,'UniformOutput',false);


tsEk = irf.ts_scalar(times,cat(2,Ek{:,iK})');

ts_Tscal = irf.ts_scalar(times,[tsT{1}.trace.data/3,tsT{2}.trace.data/3]);
ts_T_std = irf.ts_scalar(times,[std(tsT{1}.diag.data,0,2),std(tsT{2}.diag.data,0,2)]);
irf_plot({PD_use.omni.deflux.specrec,ts_Tscal,tsEk*0.008})
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
h(1).YScale = 'log';

%%
h = irf_plot(3);
Eres = 0.15;
Eratio = 1/24;

hca = irf_panel('omni');
irf_spectrogram(hca,PD.omni.deflux.specrec)

hca = irf_panel('T');
hca.ColorOrder = mms_colors('12');
irf_plot(hca,{ts_Tscal})
hca.ColorOrder = mms_colors('12');
hca.NextPlot = "add";
irf_plot(hca,{tsEk*Eres*Eratio},'--')
hca.NextPlot = "replaceall";
hca.YScale = 'log';
hca.YLabel.String = 'E (eV)';
hca.YLabel.Interpreter = "tex";
legs1 = arrayfun(@(x)sprintf('T'),1:iK,'UniformOutput',false);
legs2 = arrayfun(@(x)sprintf('%gE_k/%g',Eres,1/Eratio),1:iK,'UniformOutput',false);

hca.ColorOrder = mms_colors('1212');
legend(hca,{legs1{:},legs2{:}},'location','eastoutside')


hca = irf_panel('std(T)');
hca.ColorOrder = mms_colors('12');
irf_plot(hca,{ts_T_std})
hca.ColorOrder = mms_colors('12');
hca.YScale = 'log';
hca.YLabel.String = 'E (eV)';
legs1 = arrayfun(@(x)sprintf('std(T.diag)'),1:iK,'UniformOutput',false);
legend(hca,{legs1{:}},'location','eastoutside')

irf_plot_axis_align(h)
irf_zoom(h,'y')
h(end).XTickLabelRotation = 0;

%%
eE = logspace(1,log10(32e3),2);
%eE = PD.depend{1}(1,:)
E1 = eE(1:end-1);
E2 = eE(2:end);
dE = diff(eE);
cE = eE(1:end-1) + dE;

eV = sqrt(2*eE);
cV = sqrt(2*cE);
V1 = sqrt(2*E1);
V2 = sqrt(2*E2);
dV = V2-V1;


dEcE = dE./cE;
dVcV = dV./cV;

plot(dVcV./dEcE)

%%
E0 = 1000;
Eres = 0.14;
dE = E0*Eres;
E1 = E0-0.5*dE;
E2 = E0+0.5*dE;

v0 = sqrt(E0*2);
v1 = sqrt(E1*2);
v2 = sqrt(E2*2);
dv = v2-v1;
vres = dv/v0;

%v0 = 7;
%vres = 0.07;
N = 1000000;
dv = v0*vres;
%v = (v0-0.5*dv)+dv*vres*rand(N,1);
p = 1;
rng(1)
r = rand(N,1);
v = (v1.^p + r*(v2.^p-v1.^p)).^(1/p);
vvar = var(v); % sum((v-vmean)^2), T propto var
vmean = mean(v);
%varv/dv % the variance/v0^2 is the same for all v0^2
% T ~ 0.5*vvar ((v-vmean)^2/2)
0.5*vvar/dE

histogram(v,linspace(v1,v2,10))
%%
mu = 1000;
T_eV = 100;
sigma = T_eV*units.eV/units.mp; % m^2/s^2

dV = [dvx(:), dvy(:), dvz(:)];

% Exponent
y = dV / L;
q = sum(y.^2, 2); % L2 norm

% Maxwellian
f = n / ((2*pi)^(3/2) * prod(diag(L))) .* exp(-0.5*q); 
f = reshape(f, size(dvx));