% Generate a random initial velocity. v*f(v)

pdf = @(x) x.*exp(-x.^2);
cdf =  @(x) 1-1*exp(-x.^2);
x = linspace(0,4,100);

cdfy = rand(1e5,1);
cdfx = @(cdfy) sqrt(-log(1-cdfy));
%pdfy = pdf(cdfx);

for k = 1:4; h(k) = subplot(2,2,k); end
isub = 1;

hca=h(isub); isub=isub+1;
plot(hca,x,pdf(x),x,cdf(x))
legend(hca,'pdf','cdf')

hca=h(isub); isub=isub+1;
hist(hca,cdfy,100)

hca=h(isub); isub=isub+1;
hist(hca,cdfx(cdfy),100); hold(hca,'on');
[bincounts,binpositions] = hist(hca,cdfx(cdfy),100); hold(hca,'on');
binwidth = binpositions(2) - binpositions(1);
histarea = binwidth*sum(bincounts);

plot(hca,x,pdf(x)*histarea*2,'r'); hold(hca,'off');

hca=h(isub); isub=isub+1;
hist(hca,pdf(cdfx(cdfy)),100)


%%
% Generate a random initial velocity. v*f(v)

%R = randn(1e6,1);
nn =1e5;
vt = cn_eV2v(1600,'ev');
vd = 0;%1e5;cn_eV2v(1600,'ev');
r = vd + vt.*randn(nn,1)/sqrt(2*pi/3);
r2 = vd + vt.*randn(nn,1)/sqrt(2*pi/3);
r3 = vd + vt.*randn(nn,1);

rt = sqrt(r.^2+r2.^2);

hist(rt,100)
xlabel('v')
ylabel('pdf')
hold on

f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
%vt = 50e3; % km/s
v = linspace(-4*vt,4*vt,100);
plot(v,7e3*f(v,vt)); hold off


%%
% Generate a random initial velocity.

R = randn(1e6,1);
%rr = linspace(-3,3,100);
hist(R,100)
%%
nn = 1e4;
mu = [1 2];
Sigma = [1 .5; .5 2]; R = chol(Sigma);
z = repmat(mu,nn,1) + randn(nn,2)*R;
          
hist(z,100)

%%
nn =1e5;
vt = cn_eV2v(1600,'ev');
vd = 0;%1e5;cn_eV2v(1600,'ev');
r = vd + vt.*randn(nn,1);
r2 = vd + vt.*randn(nn,1);
r3 = vd + vt.*randn(nn,1);

rt = sqrt(r.^2+r2.^2+r3.^2);

hist(rt,100)
xlabel('v')
ylabel('pdf')
hold on
f = @(v,vt) (4/sqrt(pi))*(v.^2/(vt*sqrt(2))^3).*exp(-v.^2/(vt*sqrt(2))^2);
v = linspace(0,4*vt,100);
plot(v,1.2e3*nn*f(v,vt))
title('f = @(v,vt) (4/sqrt(pi))*(v.^2/(vt*sqrt(2))^3).*exp(-v.^2/(vt*sqrt(2))^2);')
hold off
%%
nn =1e5;
vt = cn_eV2v(1600,'ev');
vd = 0;%1e5;cn_eV2v(1600,'ev');
r = vd + vt.*randn(nn,1);
r2 = vd + vt.*randn(nn,1);
r3 = vd + vt.*randn(nn,1);

rt = sqrt(r.^2+r2.^2+0*r3.^2);

hist(rt,100)
xlabel('v')
ylabel('pdf')
hold on
f = @(v,vt) (4/sqrt(pi))*(v.^2/(vt*sqrt(2))^3).*exp(-v.^2/(vt*sqrt(2))^2);
v = linspace(0,4*vt,100);
plot(v,1.2e3*nn*f(v,vt))
title('f = @(v,vt) (4/sqrt(pi))*(v.^2/(vt*sqrt(2))^3).*exp(-v.^2/(vt*sqrt(2))^2);')
hold off
%%
vt = cn_eV2v(1600,'ev');
vd = 0;%1e5;cn_eV2v(1600,'ev');
rr = random('norm',vt,vd,[100 1] );
hist(rr,100)
xlabel('v')
ylabel('pdf')

%%
vt = cn_eV2v(1600,'ev');
f = @(v,vt) (4/sqrt(pi))*(v.^2/vt^3).*exp(-v.^2/vt^2);
v = linspace(0,4*vt,100);
plot(v,f(v,vt),v,f(v,1.5*vt),v,f(v,2*vt),v,f(v,2.5*vt))