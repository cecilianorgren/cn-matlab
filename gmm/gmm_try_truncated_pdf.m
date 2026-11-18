figure(33)
rng("default") % For reproducibility
mu0 = 0;
sigma0 = 2;
X = [];
X = [X; normrnd(mu0,sigma0,[10000,1])];
%data = [data; normrnd(mu0,sigma0,[1000,1])];
vecData = linspace(min(X),max(X),100);

xlim = [-16 16];
maxVal = 1;
heaviside_up = @(x) 1.0*(x<=0);
custompdf = @(X, mu, sigma) (normpdf(X , mu, sigma)./normcdf(maxVal, mu, sigma) .* heaviside_up(sqrt(sum(X.^2,2)) - maxVal));
%custompdf = @(X, mu, sigma) (normpdf(X , mu, sigma) .* heaviside_up(sqrt(sum(X.^2,2)) - maxVal));
customcdf = @(vmax,mu,sigma) ...
  ncx2cdf(vmax.^2./sigma.^2,3,sum(mu.^2)./sigma.^2);

keep = sqrt(sum(X.^2,2)) < maxVal;
params1 = mle(X(keep),pdf=custompdf, Start=[mu0 sigma0]);
params2 = mle(X(keep),pdf=@normpdf, Start=[mu0 sigma0]);

all.params1 = [];
all.params2 = [];
mVals = 2:14;
for mVal = mVals
  custompdf = @(X, mu, sigma) (normpdf(X , mu, sigma)./customcdf(mVal, mu, sigma) .* heaviside_up(sqrt(sum(X.^2,2)) - mVal));  

  all.params1 = [all.params1; mle(X(sqrt(sum(X.^2,2)) < mVal),pdf=custompdf, Start=[mu0 sigma0])];
  all.params2 = [all.params2; mle(X(sqrt(sum(X.^2,2)) < mVal),pdf=@normpdf, Start=[mu0 sigma0])];  
end

%
hca = subplot(2,2,1);
histogram(hca,X,Normalization="pdf")
hold(hca,'on')
plot(hca,vecData,normpdf(vecData(:),mu0,sigma0),'-')
%plot(hca,vecData,custompdf(vecData(:),params1(1),params1(2)),'-')
%plot(hca,vecData,custompdf(vecData(:),params2(1),params2(2)),'-')
hold(hca,'off')
hca.XLim = xlim;

hca = subplot(2,2,2);
histogram(hca,X(~mask),Normalization="pdf")
hold(hca,'on')
%plot(hca,vecData,normpdf(vecData(:),mu0,sigma0),'-')
plot(hca,vecData,custompdf(vecData(:),params1(1),params1(2)),'-')
plot(hca,vecData,custompdf(vecData(:),params2(1),params2(2)),'-')
plot(hca,maxVal*[1 1],hca.YLim,'--k')
hold(hca,'off')
irf_legend(hca,{"mu = " + params1(1) + ", Sigma = " + params1(2),"mu = " + params2(1) + ", Sigma = " + params2(2)}',[0.98 0.98])
hca.XLim = xlim;


%histogram(hca,data,Normalization="pdf")
%hold(hca,'on')
%plot(hca,vecData,custompdf(vecData(:),mu0,sigma0),'-')
%plot(hca,vecData,custompdf(vecData,params(1),params(2)),'-')
%hold(hca,'off')

hca = subplot(2,2,3);
plot(hca,mVals,all.params1(:,1),mVals,all.params2(:,1))
hca.XLabel.String = 'maxVal';
hca.YLabel.String = 'mu';

hca = subplot(2,2,4);
plot(hca,mVals,all.params1(:,2),mVals,all.params2(:,2))
hca.XLabel.String = 'maxVal';
hca.YLabel.String = 'Sigma';

hl = findobj(gcf,'type','line'); hl = hl(end:-1:1);
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))

