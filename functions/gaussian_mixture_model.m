function varargout = gaussian_mixture_model(R, nComp, varargin)

xvec = linspace(min(R(:,1)),max(R(:,1)),101);
yvec = linspace(min(R(:,2)),max(R(:,2)),102);
zvec = linspace(min(R(:,3)),max(R(:,3)),103);
[X,Y,Z] = ndgrid(xvec,yvec,zvec);



%initial_guess = [0 -1000 -1000;...
%                 0 -1000 +1000;...
%                 0 +1000 -1000;...
%                 0 +1000 +1000];
%S.mu = initial_guess;
%S.sigma = [500 500 500]';
%S = [];
               
%gm = fitgmdist(R,nComp,'Start',S);
gm = fitgmdist(R,nComp,varargin{:});
XYZ = [X(:) Y(:) Z(:)];
mu = gm.mu;
Sigma = gm.Sigma;
%p = mvnpdf(XYZ, mu(1,:), Sigma(:,:,1)); p = reshape(p,size(X));

gmPDF = @(x,y,z) arrayfun(@(x0,y0,z0) pdf(gm,[x0 y0 z0]),x,y,z);

for iComp = 1:nComp
  Ftmp = gm.ComponentProportion(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
  F{iComp} = reshape(Ftmp,size(X));  
end

clusterR = cluster(gm,R);

if nargout == 1
  out.R = R;
  out.idx = clusterR;
  out.gm = gm;
  out.X = X;
  out.Y = Y;
  out.Z = Z;
  out.gmPDF = gmPDF;
  out.F = F;
  varargout{1} = out;
else
  varargout{1} = clusterR;
  varargout{2} = gm;
end