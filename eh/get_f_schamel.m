function [Fsave,Ffreesave,Ftrapsave,beta] = get_f_schamel(V,n,vt,vd,PHI,VPH,nt,beta)
% [F,Ffree,Ftra,beta] = GET_F_SCHAMEL(V,n,vt,vd,PHI,VPH,nt,beta)
% if beta is a range, find bet beta witin that range, if beta is scalar,
% use that beta

units = irf_units;
E = units.me*(V-VPH).^2/2 - units.e*PHI;
dv = V(1,2) - V(1,1);

if numel(beta) == 1
  beta_vec = beta;
  n_iter = 1;
elseif numel(beta) == 2
  beta_vec = linspace(beta(1),beta(2),20);
  n_iter = 3;
else
  beta_vec = beta;
  n_iter = 1;
end

nbeta = numel(beta_vec);
C = zeros(nbeta,1);
nt_diff = Inf;
loop_iter = 0;
ibeta_best = 1;

while loop_iter < n_iter
  loop_iter = loop_iter + 1;
  for ibeta = 1:nbeta
    F = zeros(size(V));
    Ffree = zeros(size(V));
    Ftrap = zeros(size(V));
    nPop = numel(n);
    for iPop = 1:nPop
      [Ftmp,Ftmp_all] = fe_schamel(V,n(iPop),vt(iPop),vd(iPop),PHI,VPH,beta_vec(ibeta)); % distribution function, from Schamel 1986
      F = F + Ftmp;         
      Ffree = Ffree + Ftmp_all.free;
      Ftrap = Ftrap + Ftmp_all.trap;    
    end
    nt_tmp = nansum(Ftrap,2)*dv;    
    % if phi = 0, ntrap = 0, and nt_tmp will have NaN values
    nt_diff_tmp = nansum(abs(tocolumn(nt_tmp)-tocolumn(nt)));
    if 0 % plot
      hca = subplot(2,1,1);
      imagesc(hca,F')
      hold(hca,'on')
      contour(hca,E)
      hold(hca,'off')
      hcb = colorbar('peer',hca);
      hca = subplot(2,1,2);
      plot(hca,1:numel(nt),nt,1:numel(nt),nt_tmp)
      hca.Title.String = sprintf('beta = %g,  ntdiff = %g',beta_vec(ibeta),nt_diff_tmp);
      1;
    end
    if nt_diff_tmp < nt_diff
      beta = beta_vec(ibeta);
      nt_diff = nt_diff_tmp;
      Fsave = F;
      Ffreesave = Ffree;
      Ftrapsave = Ftrap;
      ibeta_best = ibeta;
    end      
  end
  if n_iter > 1
    if ibeta_best == 1
      beta_vec = linspace(beta_vec(1)-1,beta_vec(2),20);
    elseif ibeta_best == nbeta
      beta_vec = linspace(beta_vec(nbeta-1)-1,-0.001,20);
    else
      beta_vec = linspace(beta_vec(ibeta_best-1),beta_vec(ibeta_best+1),20);
    end    
    nbeta = numel(beta_vec);
  end
end
end
