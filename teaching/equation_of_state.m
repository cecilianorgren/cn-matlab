gamma = @(d) (2+d)./d;
rho = @(p,d) p.^(-gamma(d));

d = 1:3;

p = logspace(-2,0.4,100);

hca = subplot(1,1,1);
semilogy(hca,p,rho(p,d(1)),p,rho(p,d(2)),p,rho(p,d(3)));
hca.XLim =[p(1) p(end)];

hca.XLabel.String = 'p/p_0';
hca.YLabel.String = '\rho/\rho_0';
hleg = legend(hca,{sprintf('d=%g',d(1)),'d=2','d=3'});
hleg.Title.String = 'd';
hca.Title.String = '\gamma = (d+2)/2';

