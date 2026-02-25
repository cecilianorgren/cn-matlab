k = linspace(0.02,0.25,100);
lambda = 2*pi./k;
fun_lambda = @(k) 2*pi./k;

k1 = 0.05;
l1 = fun_lambda(k1);
k2 = 0.01;
l2 = fun_lambda(k2);

colors = pic_colors('matlab');

hca = subplot(1,1,1);

plot(hca,k,lambda,'displayname','\lambda = 2\pi/k')
hold(hca,'on')
plot(hca,[0 k1 k1],[l1 l1 0],'-o','displayname','k @ max theoretical growth')
plot(hca,[0 k2 k2],[l2 l2 0],'-o','displayname','k @ max observed power')
hold(hca,'off')

hca.XLabel.String = 'k';
hca.YLabel.String = '\lambda';
hca.Title.String = '';
legend(hca)
grid(hca,'on')
