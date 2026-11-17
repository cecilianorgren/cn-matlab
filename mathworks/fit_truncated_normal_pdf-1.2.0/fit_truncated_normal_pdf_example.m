%% How to fit a truncated normal (Gaussian) distribution
%%Alexey Ryabov. 2017/08/08
%% Description
%the problem is that Matlab does not allow to say that a normal distribution 
%is truncated as a result the fitted distribution will be shifted with respect
%to the data. Consider an example
clearvars;
x_min = 3;
%% Generate data normally distributed random values truncated at x_min = 1
mu1 = 5; sigma1 = 3;
pd = makedist('Normal', 'mu', mu1,'sigma', sigma1);
dat_normal = pd.random(10000, 1);
switch 'lr'
    case 'l'
        dat_normal = dat_normal(dat_normal>x_min);
    case 'r'
        dat_normal = dat_normal(dat_normal<10);
    case 'lr'
        dat_normal = dat_normal(dat_normal>x_min & dat_normal<10);
end

%% Make a historgam
x = -5:0.2:25;
[counts,centers]  = hist(dat_normal, x);
bar(centers, counts/(sum(counts) * mean(diff(x)) ))
hold on
%use fitdist to find the paramters of our distribution
pd_fit = fitdist(dat_normal,'Normal');
display(pd_fit)

[norm_trunc, phat, phat_ci] = fitdist_ntrunc(dat_normal);
%% Plot results
figure(1)
clf
%Data
bar(centers, counts/(sum(counts) * mean(diff(x)) ), 0.4, 'Edgecolor', 'none');
hold on

%True distribution
plot(x, (norm_trunc(x , mu1 , sigma1)), 'Color', [0.47, 0.67, 0.19], 'Linewidth', 5)

%the distribution obtained with fitdist is shifted to the right
plot(x, pd_fit.pdf(x),  'Linewidth', 3);

%fitted truncated distribution
plot(x, (norm_trunc(x , phat(1), phat(2))), '--', 'Color', [237, 175, 33]/255, 'Linewidth', 3)
hold off
legend({'Data', 'True distribution', 'fit Gaussian' 'fit truncated'})
xlabel('x');
ylabel('Probability density')

