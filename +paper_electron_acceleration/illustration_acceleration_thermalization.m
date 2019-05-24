f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt.^2); % this is reduced

h = setup_subplots(1,1);
isub = 1;
hca = h(isub); isub = isub + 1;

linewidth = 1.5;

vmin = -20000*1e3;
vmax = 50000*1e3;
nv = 1000;
vvec = linspace(vmin,vmax,nv);

[n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('ill_lobe');
nsp = numel(n);
fplot = zeros(1,nv);
for isp = 1:nsp, fplot = fplot + f(vvec,n(isp),vt(isp),vd(isp)); end
plot(hca,vvec*1e-6,fplot,'linewidth',linewidth,'displayname','lobe')

hold(hca,'on')

[n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('ill_lobe_acc');
nsp = numel(n);
fplot = zeros(1,nv);
for isp = 1:nsp, fplot = fplot + f(vvec,n(isp),vt(isp),vd(isp)); end
plot(hca,vvec*1e-6,fplot,'linewidth',linewidth,'displayname','accelerated lobe - beam')

[n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('ill_lobe_acc_therm1');
nsp = numel(n);
fplot = zeros(1,nv);
for isp = 1:nsp, fplot = fplot + f(vvec,n(isp),vt(isp),vd(isp)); end
plot(hca,vvec*1e-6,fplot,'linewidth',linewidth,'displayname','partly thermalized beam')

[n,T,m,q,vd,vt,wp,Lin,Ld] = f_inp('ill_lobe_acc_therm2');
nsp = numel(n);
fplot = zeros(1,nv);
for isp = 1:nsp, fplot = fplot + f(vvec,n(isp),vt(isp),vd(isp)); end
plot(hca,vvec*1e-6,fplot,'linewidth',linewidth,'displayname','partly thermalized beam')

hold(hca,'off')

hca.XLabel.String = 'v';
hca.YLabel.String = 'f_e';
hca.XTick = [];
hca.YTick = [];
hca.XLim = [vmin vmax]*1e-6;
hca.FontSize = 12;

hleg = legend(hca,'location','northeast');
hleg.Box = 'off';