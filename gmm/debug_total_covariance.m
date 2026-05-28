w = [0.5 0.5];
mu = [00 0 0; 700 0 0];
Sig = repmat(eye(3,3),[1 1 2])*200^2;
[wS,muS,SigS] = total_covariance(w,mu,Sig);

vvec = -2000:50:2000;
[X,Y,Z] = ndgrid(vvec,vvec,vvec);

XYZ = [X(:) Y(:) Z(:)];

Ftot = zeros(size(X));
Fcomp = zeros([size(X) numel(w)]);

Fmerg = wS*mvnpdf(XYZ, muS, SigS);
Fmerg = reshape(Fmerg,size(X));  

for iComp = 1:numel(w)
  Ftmp = w(iComp)*mvnpdf(XYZ, mu(iComp,:), Sig(:,:,iComp));
  Ftmp = reshape(Ftmp,size(X));  
  Ftot = Ftot + Ftmp;
  Fcomp(:,:,:,iComp)  = Ftmp; 
end

% Sample particles from distributions
g_o.ComponentProportion = w;
g_o.mu = mu;
g_o.Sigma = Sig;

g_m.ComponentProportion = wS;
g_m.mu = muS;
g_m.Sigma = SigS;

N = 100000;
xyz_p_o =  gmm_monte_carlo_sampling(g_o,N,1:2);
xyz_p_o1 = gmm_monte_carlo_sampling(g_o,round(N*w(1)/sum(w)),1);
xyz_p_o2 = gmm_monte_carlo_sampling(g_o,round(N*w(2)/sum(w)),2);
xyz_p_m =   gmm_monte_carlo_sampling(g_m,N,1);

idim = 1;
[pdf_p_o,~,mid] = histcn(xyz_p_o(:,idim), vvec);
[pdf_p_o1,~,mid] = histcn(xyz_p_o1(:,idim), vvec);
[pdf_p_o2,~,mid] = histcn(xyz_p_o2(:,idim), vvec);
[pdf_p_m,~,mid] = histcn(xyz_p_m(:,idim), vvec);

% Plot results
hca = subplot(2,1,1);
f_merg = squeeze(sum(Fmerg,[2 3]));
f_comp = squeeze(sum(Fcomp,[2 3]));
f_tot = squeeze(sum(Ftot,[2 3]));
plot(hca, vvec,f_comp, vvec,f_tot, vvec,f_merg)

hca = subplot(2,1,2);
plot(hca, mid{1},pdf_p_o1, mid{1},pdf_p_o2, mid{1},pdf_p_o, mid{1},pdf_p_m)

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
linkprop(h,{'XLim'});