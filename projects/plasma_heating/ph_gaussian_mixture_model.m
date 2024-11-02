%% 2D
mu1 = [1 2];          % Mean of the 1st component
sigma1 = [2 0; 0 .5]; % Covariance of the 1st component
mu2 = [-3 -5];        % Mean of the 2nd component
sigma2 = [1 0; 0 1];  % Covariance of the 2nd component

rng('default') % For reproducibility
r1 = mvnrnd(mu1,sigma1,1000);
r2 = mvnrnd(mu2,sigma2,1000);
X = [r1; r2];

xvec = linspace(min(R(:,1)),max(R(:,1)),101);
yvec = linspace(min(R(:,2)),max(R(:,2)),102);
[Xv,Yv] = meshgrid(xvec,yvec);

gm = fitgmdist(X,2);

XXv = [Xv(:) Yv(:)];
mu = squeeze(gm.mu(1,:));
Sigma = squeeze(gm.Sigma(:,:,1));
p = mvnpdf(XXv, mu, Sigma);
p = reshape(p,numel(xvec),numel(yvec));

scatter(X(:,1),X(:,2),10,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-8 6])
%% 2D
mu1 = [1 2];          % Mean of the 1st component
sigma1 = [2 0; 0 .5]; % Covariance of the 1st component
mu2 = [-3 -5];        % Mean of the 2nd component
sigma2 = [1 0; 0 1];  % Covariance of the 2nd component

rng('default') % For reproducibility
r1 = mvnrnd(mu1,sigma1,1000);
r2 = mvnrnd(mu2,sigma2,1000);
X = [r1; r2];

xvec = linspace(min(R(:,1)),max(R(:,1)),101);
yvec = linspace(min(R(:,2)),max(R(:,2)),102);
[Xv,Yv] = meshgrid(xvec,yvec);

gm = fitgmdist(X,2);

XXv = [Xv(:) Yv(:)];
mu = squeeze(gm.mu(1,:));
Sigma = squeeze(gm.Sigma(:,:,1));
p = mvnpdf(XXv, mu, Sigma);
p = reshape(p,numel(xvec),numel(yvec));

scatter(X(:,1),X(:,2),10,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-8 6])
%% 3D
mu1 = [1 2 0];          % Mean of the 1st component
sigma1 = [2 0 0; 0 .5 00; 0 0 0.5]; % Covariance of the 1st component
mu2 = [-3 -5 1];        % Mean of the 2nd component
sigma2 = [1 0 0; 0 1 0; 0 0 0.1];  % Covariance of the 2nd component
mu3 = [-3 2 -1];        % Mean of the 3rd component
sigma3 = [1 0 0; 0 1 0; 0 0 0.2];  % Covariance of the 3rd component



rng('default') % For reproducibility
r1 = mvnrnd(mu1,sigma1,1000);
r2 = mvnrnd(mu2,sigma2,1000);
r3 = mvnrnd(mu3,sigma3,1000);
R = [r1; r2; r3];
%%
%R = V;

xvec = linspace(min(R(:,1)),max(R(:,1)),101);
yvec = linspace(min(R(:,2)),max(R(:,2)),102);
zvec = linspace(min(R(:,3)),max(R(:,3)),103);
[X,Y,Z] = ndgrid(xvec,yvec,zvec);

nComp = 4;


initial_guess = [0 -1000 -1000;...
                 0 -1000 +1000;...
                 0 +1000 -1000;...
                 0 +1000 +1000];
S.mu = initial_guess;
S.sigma = [500 500 500]';
S = [];
               
%gm = fitgmdist(R,nComp,'Start',S);
gm = fitgmdist(R,nComp);
XYZ = [X(:) Y(:) Z(:)];
mu = gm.mu;
Sigma = gm.Sigma;
%p = mvnpdf(XYZ, mu(1,:), Sigma(:,:,1)); p = reshape(p,size(X));

gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm,[x0 y0 z0]),x,y,z);
F1 = gm.ComponentProportion(1)*mvnpdf(XYZ, mu(1,:), Sigma(:,:,1)); F1 = reshape(F1,size(X));
F2 = gm.ComponentProportion(2)*mvnpdf(XYZ, mu(2,:), Sigma(:,:,2)); F2 = reshape(F2,size(X));
F3 = gm.ComponentProportion(3)*mvnpdf(XYZ, mu(3,:), Sigma(:,:,3)); F3 = reshape(F3,size(X));

clusterR = cluster(gm,R);
%%
hca = subplot(1,1,1);

%scatter3(R(:,1),R(:,2),R(:,3),20,'.') % Scatter plot with points of size 10
scatter3(R(:,1),R(:,2),R(:,3),20,clusterR) % Scatter plot with points of size 10
hold(hca,'on')


F = gmPDF(X,Y,Z);


% Flev = iso_values(isurf);
s = isosurface(X,Y,Z,F,2*1e-3);
s1 = isosurface(X,Y,Z,F1,2*1e-3);
s2 = isosurface(X,Y,Z,F2,2*1e-3);
s3 = isosurface(X,Y,Z,F3,2*1e-3);

%p = patch(hca,'Faces',s.faces,'Vertices',s.vertices,'FaceAlpha',0.5);
p1 = patch(hca,'Faces',s1.faces,'Vertices',s1.vertices,'FaceAlpha',0.2,'FaceColor','r','EdgeColor','none');
p2 = patch(hca,'Faces',s2.faces,'Vertices',s2.vertices,'FaceAlpha',0.2,'FaceColor','b','EdgeColor','none');
p3 = patch(hca,'Faces',s3.faces,'Vertices',s3.vertices,'FaceAlpha',0.2,'FaceColor','g','EdgeColor','none');

lighting(hca,'gouraud')
camlight(hca);

hold(hca,'off')

%%
x = [randn(4000,1)/2; 5+2*randn(6000,1)];
f = fitgmdist(x,2);
histogram(x,'Normalization','pdf')
xgrid = linspace(-4,12,1001)';
hold on; plot(xgrid,pdf(f,xgrid),'r-'); hold off

n1 = makedist('normal',f.mu(1),sqrt(f.Sigma(1)));
n2 = makedist('normal',f.mu(2),sqrt(f.Sigma(2)));
p = f.ComponentProportion;
y1 = p(1)*pdf(n1,xgrid);
y2 = p(2)*pdf(n2,xgrid)
hold on; plot(xgrid,y1,xgrid,y2,'c--'); hold off

%%
load fisheriris;
X = meas(:,1:2);
[n,p] = size(X);

plot(X(:,1),X(:,2),'.','MarkerSize',15);
title('Fisher''s Iris Data Set');
xlabel('Sepal length (cm)');
ylabel('Sepal width (cm)');

rng(3);
k = 3; % Number of GMM components
options = statset('MaxIter',1000);

Sigma = {'diagonal','full'}; % Options for covariance matrix type
nSigma = numel(Sigma);

SharedCovariance = {true,false}; % Indicator for identical or nonidentical covariance matrices
SCtext = {'true','false'};
nSC = numel(SharedCovariance);

d = 500; % Grid length
x1 = linspace(min(X(:,1))-2, max(X(:,1))+2, d);
x2 = linspace(min(X(:,2))-2, max(X(:,2))+2, d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];

threshold = sqrt(chi2inv(0.99,2));
count = 1;
for i = 1:nSigma
    for j = 1:nSC
        gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i}, ...
            'SharedCovariance',SharedCovariance{j},'Options',options); % Fitted GMM
        clusterX = cluster(gmfit,X); % Cluster index 
        mahalDist = mahal(gmfit,X0); % Distance from each grid point to each GMM component
        % Draw ellipsoids over each GMM component and show clustering result.
        subplot(2,2,count);
        h1 = gscatter(X(:,1),X(:,2),clusterX);
        hold on
            for m = 1:k
                idx = mahalDist(:,m)<=threshold;
                Color = h1(m).Color*0.75 - 0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
                uistack(h2,'bottom');
            end    
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s\nSharedCovariance = %s',Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'1','2','3'})
        hold off
        count = count + 1;
    end
end


%% mvnpdf
mu = [1 -1]; Sigma = [.9 .4; .4 .3];
[X1,X2] = meshgrid(linspace(-1,3,25)', linspace(-3,1,25)');
X = [X1(:) X2(:)];
p = mvnpdf(X, mu, Sigma);
surf(X1,X2,reshape(p,25,25));

%% Test on MMS distribution


 