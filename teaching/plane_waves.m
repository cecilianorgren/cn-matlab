wpe = 1;

A = @(k,x,t) cos(k.*x-wpe*t)+cos(k.*x+wpe*t);

k = 1;
lambda = 2*pi/k;
nt = 10;
t = wpe^-1*linspace(0,4,nt);
nx = 100;
x = lambda*linspace(0,2,nx);

hca = subplot(1,1,1);
for it = 1:nt  
  plot(hca,x,A(k,x,t(it)))
  hca.YLim = 2*[-1 1];
  pause(0.1)
end


%% Ionospheric reflection
c = 1;
w0 = 1; % this is our launched wave
k0 = w0/c;

k = @(w0,wpe) sqrt(w0.^2-wpe.^2)/c; % wavenumber of our launched wave
w_ = @(k,wpe) sqrt(wpe.^2 + c^2.*k.^2); % general dispersion relation

nwpe = 7; % error if it exceeds number of colors
wpe = w0*linspace(0,1,nwpe);

nk = 100;
k_ = linspace(0,2,nk);

hca = subplot(1,1,1);
colors = get(hca,'colororder');
plot(hca,0,0) % just so the hold on works
hold(hca,'on')
clear hleg legends
for iwpe = 1:nwpe
  hleg_w(iwpe) = plot(hca,k_,w_(k_,wpe(iwpe)),'color',colors(iwpe,:));
  legends{iwpe,1} = sprintf('w_{pe}/w_0 = %.2f',wpe(iwpe));
  hleg_w0(iwpe) = plot(hca,k(w0,wpe(iwpe)),w_(k(w0,wpe(iwpe)),wpe(iwpe)),'color',colors(iwpe,:),'marker','o');
  hleg_wpe(iwpe) = plot(hca,k_,k_*0+wpe(iwpe),'--','color',colors(iwpe,:));
end
hold(hca,'off')
hca.XLabel.String = 'k/k_0';
hca.YLabel.String = '\omega/\omega_0';
legend(hca,[hleg_w(1),hleg_w0(1),hleg_wpe(1)],{'\omega(k)','\omega_0','\omega_{pe}'},'location','north')
%a=axes('position',get(hca,'position'),'visible','off');
%legend(a,[hleg_w],legends,'location','northwest')
irf_legend(hca,legends,[0.01 0.98])

