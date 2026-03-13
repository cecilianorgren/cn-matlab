
a = randn(30,1)/2;
b = randn(50,1)/2 + 3;
data = [a; b];
n = numel(data);

K = 3;
gm = fitgmdist(data,K);

likelihood_comp = zeros(n,K);
for ic = 1:K
  %for ip = 1:n
    p{ic} = mvnpdf(data, gm.mu(ic,:), gm.Sigma(:,:,ic));
    likelihood_comp(:,ic) = gm.ComponentProportion(ic)*p{ic};
    %nlogL{ic} = normlike([gm.mu(ic),gm.Sigma(ic)],data);
  %end
end


likelihood = sum(likelihood_comp,2); % sum over components to get p(x_i)
logL = sum(log(likelihood));


colors = pic_colors('matlab');

h = setup_subplots(2,1,'vertical');
isub = 1;

hca = h(isub); isub = isub + 1;
%[N,edges] = histcounts(data,20);
histogram(hca,data,20,'Normalization','pdf','FaceColor',[0.4 0.4 0.4])
hca.NextPlot = "add";
hca.ColorOrder = [0.4 0.4 0.4; colors];
nvec = linspace(min(data),max(data),100);
for ic = 1:K
  toplot = gm.ComponentProportion(ic)*mvnpdf(nvec',gm.mu(ic,:), gm.Sigma(:,:,ic));
  plot(hca,nvec,toplot,'linewidth',2)
end
hca.NextPlot = "replace";
hca.XLabel.String = 'Data value';
hca.YLabel.String = 'PDF(x_i)';


if 0 % histogram of likelihood
hca = h(isub); isub = isub + 1;
histogram(hca,likelihood_comp(:,1),20,'Normalization','pdf')
hca.NextPlot = "add";
histogram(hca,likelihood_comp(:,2),20,'Normalization','pdf')
hca.NextPlot = "replace";
end

%hca = h(isub); isub = isub + 1;
%plot(hca,likelihood(:,1),likelihood(:,2),'x')


hca = h(isub); isub = isub + 1;
%hca.ColorOrder = colors([3 1 2],:);
[~,idx] = sort(likelihood_comp(:,1));
%[a,p1,p2] = plotyy(hca,likelihood_comp(idx,:),'o',likelihood(idx),'x');
plot(hca,likelihood_comp(idx,:),'.');
hca.NextPlot = "add";
plot(hca,likelihood(idx),'o');
hca.NextPlot = "replace";
hca.XLabel.String = 'Data point, i';
hca.YLabel.String = 'Probability';
%a(2).YLabel.String = {'Likelihood','(Sum of probabilites)'};
hca.Title.String = 'Data sorted by ascending likelihood of component 1';
%plot(hca,likelihood(idx,:),'.')
%hca.NextPlot = "add";
%plotyy(hca,prod(likelihood(idx,:),2),'o',prod(likelihood(idx,:),2),'x')
%hca.NextPlot = "replace";
legs = arrayfun(@(x) sprintf('p_%g(x_i)',x),1:K,'UniformOutput',false);
irf_legend(hca,{legs{:},'\Sigma_{k=1:K} p_k(x_i)'}',[0.02 0.98])


if 0
hca = h(isub); isub = isub + 1;
[~,idx] = sort(likelihood(:,1));
plot(hca,log(likelihood_comp(idx,:)),'.');
hca.NextPlot = "add";
plot(hca,log(likelihood(idx)),'o')
hca.NextPlot = "replace";
hca.XLabel.String = 'Data point, i';
hca.YLabel.String = 'log(Probability)';
hca.Title.String = 'Data sorted by ascending probability of component 1';
irf_legend(hca,{sprintf('own calculated -log(likelihood) = %g',-logL),sprintf('gm.NegativeLogLikelihood = %g',gm.NegativeLogLikelihood)}',[0.98 0.02])
end

if 0
hca = h(isub); isub = isub + 1;
[~,idx] = sort(likelihood(:,1));
plot(hca,prod(likelihood(idx,:),2),'.')
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 14;',1:numel(h))

%-logL/gm.NegativeLogLikelihood % log10(e) = 0.4343, where does this difference come from, because I used log10...
gm.BIC