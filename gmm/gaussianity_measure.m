function varargout = gaussianity_measure(g_o,N,vvec)

if ~exist('N','var')
  N = 10000;
end

g_o.NumComponents = numel(g_o.w);
g_o.ComponentProportion = g_o.w;

g_m = total_covariance(g_o.w,g_o.mu,g_o.Sigma);
g_m.ComponentProportion = g_m.w;
g_m.NumComponents = 1;


%vvec = -2000:50:2000;
[X,Y,Z] = ndgrid(vvec,vvec,vvec);
XYZ = [X(:) Y(:) Z(:)];

Ftot = zeros(size(X));
Fcomp = zeros([size(X) numel(g_o.w)]);

Fmerg = g_m.w*mvnpdf(XYZ, g_m.mu, g_m.Sigma);
Fmerg = reshape(Fmerg,size(X));  

for iComp = 1:numel(g_o.w)
  Ftmp = g_o.w(iComp)*mvnpdf(XYZ, g_o.mu(iComp,:), g_o.Sigma(:,:,iComp));
  Ftmp = reshape(Ftmp,size(X));  
  Ftot = Ftot + Ftmp;
  Fcomp(:,:,:,iComp)  = Ftmp; 
end

p_grid = squeeze(sum(Fcomp(:,:,:,1:2),4));
q_grid = Fmerg;

% Sample particles from distributions
%g_o
%g_m
xyz_p_o =  gmm_monte_carlo_sampling(g_o,N,1:g_o.NumComponents);
xyz_p_m =  gmm_monte_carlo_sampling(g_m,N,1);

% Evaluate gaussianity
p = zeros(N,1);
for ik = 1:g_o.NumComponents
  ptmp = g_o.w(ik)*mvnpdf(xyz_p_o, g_o.mu(ik,:), g_o.Sigma(:,:,ik));
  p = p + ptmp;
end

q = sum(g_m.w)*mvnpdf(xyz_p_o, g_m.mu(1,:), g_m.Sigma(:,:,1));

%[~, ~, p_] = gmm_evaluate_f(g_o,xyz_p_o,1,'group',{[1 2]});
%[q_, ~] = gmm_evaluate_f(g_m,xyz_p_o,1); % evaluate f at xyz
% Calculate rms
rms_mc = sqrt(sum((p-q).^2,'all'))/sqrt(sum(p.^2,'all'));
rms_grid = sqrt(sum((p_grid-q_grid).^2,'all'))/sqrt(sum(p_grid.^2,'all'));



varargout{1} = rms_mc;
varargout{2} = rms_grid;




if 0 % plot
idim_vec = 1:3;

h = setup_subplots(numel(idim_vec),2);
isub = 1;

for idim = idim_vec
  
  [pdf_p_o,~,mid] = histcn(xyz_p_o(:,idim), vvec);
  [pdf_p_o1,~,mid] = histcn(xyz_p_o1(:,idim), vvec);
  [pdf_p_o2,~,mid] = histcn(xyz_p_o2(:,idim), vvec);
  [pdf_p_m,~,mid] = histcn(xyz_p_m(:,idim), vvec);
  
  % Plot results
  hca = h(isub); isub = isub + 1;
  f_merg = squeeze(sum(Fmerg,setdiff(1:3,idim)));
  f_comp = squeeze(sum(Fcomp,setdiff(1:3,idim)));
  f_tot = squeeze(sum(Ftot,setdiff(1:3,idim)));
  plot(hca, vvec,f_comp, vvec,f_tot, vvec,f_merg)
  hca.XLabel.String = 'v (km/s)';
  hca.Title.String = 'Grid';
  legend(hca,{'f^1','f^2','f^{1+2} (p)','f^{merge} (q)'},'box','off')
  
  hca = h(isub); isub = isub + 1;
  plot(hca, mid{1},pdf_p_o1, mid{1},pdf_p_o2, mid{1},pdf_p_o, mid{1},pdf_p_m)
  hca.XLabel.String = 'v (km/s)';
  hca.Title.String = 'Montecarlo';
  legend(hca,{'f^1','f^2','f^{1+2} (p)','f^{merge} (q)'},'box','off')
  irf_legend(hca,{'RMS = sum(sqrt((q-p)^2)/sqrt((p)^2))',sprintf('RMS montecarlo = %0.4f',rms_mc),sprintf('RMS grid = %0.4f',rms_grid)}',[0.02 0.98],'color','k')
end

%h = findobj(gcf,'type','axes'); h = h(end:-1:1);
linkprop(h,{'XLim'});
end