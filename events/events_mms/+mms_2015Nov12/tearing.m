%% Liu2013, guide field, force free magnetic field
units = irf_units;

Delta = @(k,lambda,bg,theta) 2*(1+bg^2*tand(theta).^2)./k/lambda^2 - 2*k;
ls = @(lambda,bg,theta) (cosd(theta).^2-bg^2.*sind(theta).^2)./(lambda*cosd(theta)*sqrt(1+bg^2));
wpe = @(ne) sqrt(ne*units.e^2/units.eps0/units.me);
de = @(ne) units.c/wpe(ne);
vthe = @(Te) sqrt(2*units.eV*Te/units.me);
gamma = @(k,lambda,bg,theta,Te,Ti,ne) k*vthe(Te)*de(ne)^2.*Delta(k,lambda,bg,theta)./(2*sqrt(pi)*ls(lambda,bg,theta)*(1+sqrt(units.me*Te/units.mp/Ti)));

bg = 0.0;
Te = 50;
Ti = Te*5;500;
theta = 0:90;
ne = 5e6;

nRows = 4;
nCols = 1;
for iSubPlot = 1:nRows*nCols;
  h(iSubPlot) = subplot(nRows,nCols,iSubPlot);
end

isub = 1;
bg = 0;

% gamma as a function of k*lambda
hca = h(isub); isub = isub + 1;
theta = [0 5 10];
lambda = 2*de(ne);
k = linspace(0.01,1.4,100)/lambda;

plot(hca,k*lambda,[gamma(k,lambda,bg,theta(1),Te,Ti,ne*1e6); gamma(k,lambda,bg,theta(2),Te,Ti,ne*1e6); gamma(k,lambda,bg,theta(3),Te,Ti,ne*1e6)])
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '\gamma';
legend(hca,{'\theta = 0','\theta = 5','\theta = 10'})
hca.YGrid = 'on';

% gamma as a function of theta
hca = h(isub); isub = isub + 1;
lambda = 2*de(ne);
k = 0.5/lambda; % k*lamda = 0.5;
theta = 1:2:90;
gamma0 = gamma(k,lambda,bg,0,Te,Ti,ne*1e6);
plot(hca,theta,gamma(k,lambda,bg,theta,Te,Ti,ne*1e6)/gamma0);
hca.XLabel.String = '\theta';
hca.YLabel.String = '\gamma/\gamma_0';

% Delta
hca = h(isub); isub = isub + 1;
theta = 0;
lambda = 2*de(ne);
k = linspace(0.001,1.4,100)/lambda;
%plotDelta = [Delta(k,lambda,bg,theta(1)); Delta(k,lambda,bg,theta(2)); Delta(k,lambda,bg,theta(3))];

%plotDelta(plotDelta<0) = NaN;
plot(hca,k*lambda,k.*Delta(k,lambda,bg,theta));
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '\Delta''';

% 1/ls
hca = h(isub); isub = isub + 1;
theta = 0:1:90;
lambda = 2*de(ne);
k = 0.5/lambda;
plotyy(hca,theta,ls(lambda,bg,theta),theta,1./ls(lambda,bg,theta));
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '1/ls';

%% Daughton2011 (Nature), guide field
units = irf_units;

wp = @(n,m) sqrt(n*units.e^2/units.eps0/m);
d = @(n,m) units.c/wp(n,m);
vth = @(T,m) sqrt(2*units.eV*T/m);
wc = @(B,m) units.e*B/m;
rho = @(T,m,B) vth(T,m)./wc(B,m);

Bx0 = 10*1e-9; % T
By0 = 5*1e-9;
bg = By0/Bx0;
Te = 50; % eV
Ti = 10*Te;
ne = 5e6;
nb = 0.3;

zs = @(lambda,bg,theta) -lambda*atanh(tand(theta)*bg);
alpha = @(Te,Ti) sqrt(Te*units.me/Ti/units.mp);
%Delta = @(k,lambda,bg,theta) 2/lambda*(1./k./lambda-k.*lambda).*(1+0.5*tanh(zs(lambda,theta,bg)).^2.*(1+1./(1-k*lambda)));
Delta = @(k,lambda,bg,theta) 2./lambda.*(1./k./lambda-k.*lambda).*(1+0.5*tand(theta).^2*bg.^2*(1+1./(1-k*lambda)));
G = @(k,lambda,bg,theta,nb) (1-sind(theta).^2*(1-bg.^2*(1-0.5*k*lambda)./(1-k*lambda)))*1./(1+nb*cosh(zs(lambda,bg,theta)/lambda).^2);
wi = @(k,lambda,Te,Ti,B,nb,theta) wc(B,units.mp)*sqrt(1/pi)*alpha(Te,Ti)./(1+alpha(Te,Ti)).*(rho(Ti,units.mp,B)/lambda).^3.*(1+Te./Ti).*(1-k.^2*lambda^2).*G(k,lambda,bg,theta,nb)./bg;



nRows = 4;
nCols = 1;
for iSubPlot = 1:nRows*nCols;
  h(iSubPlot) = subplot(nRows,nCols,iSubPlot);
end

isub = 1;

lambda = 5*d(ne,units.me); % initial half thickness of current layer
%lambda = 1*rho(Ti,units.mp,Bx0);
kl = 0.4;

% G as a function of theta
hca = h(isub); isub = isub + 1;
theta = 0:0.01:90;
k = kl/lambda;
plotG = [G(k,lambda,bg,theta,0.01); G(k,lambda,bg,theta,0.1); G(k,lambda,bg,theta,0.5); G(k,lambda,bg,theta,0.9)];
maxG = max(G(k,lambda,bg,theta,0.0));
plotG(plotG<0) = NaN;
plotG(plotG>maxG) = NaN;
plot(hca,theta,plotG); %hold(hca,'on')
%plot(hca,theta,[G(k,lambda,bg,theta,0.1); G(k,lambda,bg,theta,0.5); G(k,lambda,bg,theta,0.9)],'--'); hold(hca,'off') 
hca.XLabel.String = '\theta';
hca.YLabel.String = 'G';
legend(hca,{'n_b/n_0 = 0.01','n_b/n_0 = 0.1','n_b/n_0 = 0.5','n_b/n_0 = 0.9'},'location','best')
hca.YGrid = 'on';
hca.Title.String = sprintf('b_g = B_{y0}/B_{x0} = %g',bg);

% G as a function of k*lambda
hca = h(isub); isub = isub + 1;
theta = [10 20];
k = linspace(0.001,0.9,100)/lambda;
plotG = [G(k,lambda,bg,theta(1),nb);G(k,lambda,bg,theta(2),nb)];
plotG(plotG<0) = NaN;
plot(hca,k*lambda,plotG);
legend(hca,{['\theta = ' num2str(theta(1))],['\theta = ' num2str(theta(2))]},'location','best')
hca.XLabel.String = 'k*\lambda';
hca.YLabel.String = 'G';

% gamma as a function of k*lambda
hca = h(isub); isub = isub + 1;
theta = [10 20];
k = linspace(0.001,1.4,100)/lambda;
plotWi = [wi(k,lambda,Te,Ti,Bx0,nb,theta(1));wi(k,lambda,Te,Ti,Bx0,nb,theta(2))]/wc(Bx0,units.mp);
plot(hca,k*lambda,plotWi);
legend(hca,{['\theta = ' num2str(theta(1))],['\theta = ' num2str(theta(2))]})
hca.XLabel.String = 'k*\lambda';
hca.YLabel.String = '\gamma/\omega_{ci}';

% gamma as a function of theta, for fixed k*lambda
hca = h(isub); isub = isub + 1;
k = kl/lambda;
theta = 1:0.5:90;
lambda_rel = [1 2]*lambda;

plotWi = [wi(k,lambda_rel(1),Te,Ti,Bx0,nb,theta); wi(k,lambda_rel(2),Te,Ti,Bx0,nb,theta)]/wc(Bx0,units.mp);
for ii = 1:size(plotWi,1)
  tmp_max = findpeaks(real(plotWi(ii,:)));
  maxWi(ii) = tmp_max(1);
end
minWi = min(min(abs(real(plotWi))));

semilogy(hca,theta,plotWi);
hca.XLabel.String = '\theta';
hca.YLabel.String = '\gamma/\omega_{ci}';
legend(hca,{['\lambda =  ' num2str(lambda_rel(1)/lambda) '\lambda_0'],['\lambda =  ' num2str(lambda_rel(2)/lambda) '\lambda_0']})
hca.YLim(1) = 10.^floor(log10(minWi));
hca.YLim(2) = 10.^ceil(log10(max(maxWi)));
hca.YTick = 10.^(log10(hca.YLim(1)):log10(hca.YLim(2)));
hca.YGrid = 'on';
hca.YMinorGrid = 'off';
hca.Title.String = sprintf('b_g = B_{y0}/B_{x0} = %g, T_i/T_e = %g, n_b/n_0 = %g, lambda = %.0f d_e = %.3g rho_{i} (B_{x0}= %d nT)',bg,Ti/Te,nb,lambda/d(ne,units.me),lambda/rho(Ti,units.mp,Bx0),Bx0*1e9);

%% Karimabadi2005b, guide field
units = irf_units;

Delta = @(k,lambda,bg,theta) 2*(1+bg^2*tand(theta).^2)./k/lambda^2 - 2*k;
ls = @(lambda,bg,theta) (cosd(theta).^2-bg^2.*sind(theta).^2)./(lambda*cosd(theta)*sqrt(1+bg^2));
wpe = @(ne) sqrt(ne*units.e^2/units.eps0/units.me);
de = @(ne) units.c/wpe(ne);
vthe = @(Te) sqrt(2*units.eV*Te/units.me);
gamma = @(k,lambda,bg,theta,Te,Ti,ne) k*vthe(Te)*de(ne)^2.*Delta(k,lambda,bg,theta)./(2*sqrt(pi)*ls(lambda,bg,theta)*(1+sqrt(units.me*Te/units.mp/Ti)));

bg = 0.0;
Te = 50;
Ti = Te*5;500;
theta = 0:90;
ne = 5e6;

nRows = 4;
nCols = 1;
for iSubPlot = 1:nRows*nCols;
  h(iSubPlot) = subplot(nRows,nCols,iSubPlot);
end

isub = 1;
bg = 0;

% gamma as a function of k*lambda
hca = h(isub); isub = isub + 1;
theta = [0 5 10];
lambda = 2*de(ne);
k = linspace(0.01,1.4,100)/lambda;

plot(hca,k*lambda,[gamma(k,lambda,bg,theta(1),Te,Ti,ne*1e6); gamma(k,lambda,bg,theta(2),Te,Ti,ne*1e6); gamma(k,lambda,bg,theta(3),Te,Ti,ne*1e6)])
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '\gamma';
legend(hca,{'\theta = 0','\theta = 5','\theta = 10'})
hca.YGrid = 'on';

% gamma as a function of theta
hca = h(isub); isub = isub + 1;
lambda = 2*de(ne);
k = 0.5/lambda; % k*lamda = 0.5;
theta = 1:2:90;
gamma0 = gamma(k,lambda,bg,0,Te,Ti,ne*1e6);
plot(hca,theta,gamma(k,lambda,bg,theta,Te,Ti,ne*1e6)/gamma0);
hca.XLabel.String = '\theta';
hca.YLabel.String = '\gamma/\gamma_0';

% Delta
hca = h(isub); isub = isub + 1;
theta = 0;
lambda = 2*de(ne);
k = linspace(0.001,1.4,100)/lambda;
%plotDelta = [Delta(k,lambda,bg,theta(1)); Delta(k,lambda,bg,theta(2)); Delta(k,lambda,bg,theta(3))];

%plotDelta(plotDelta<0) = NaN;
plot(hca,k*lambda,k.*Delta(k,lambda,bg,theta));
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '\Delta''';

% 1/ls
hca = h(isub); isub = isub + 1;
theta = 0:1:90;
lambda = 2*de(ne);
k = 0.5/lambda;
plotyy(hca,theta,ls(lambda,bg,theta),theta,1./ls(lambda,bg,theta));
hca.XLabel.String = 'k\lambda';
hca.YLabel.String = '1/ls';
