function out = merge_gaussians_maxwellianity(gm,varargin)

% Defaults
doSort = false;
method = 'montecarlo';
N = 100000;
vvec = -2500:50:2500;

nargs = numel(varargin);
if nargs > 0, have_options = 1; args = varargin(:); end

while have_options
    switch lower(args{1})
      case {'method'}
        method = args{2};               
        args = args(3:end);
      case {'nmc'}
        N = args{2};
        args = args(3:end);
      case {'grid'}
        vvec = args{2};
        args = args(3:end);
      case 'sort'
        doSort = args{2};
        args = args(3:end);
      otherwise
        error('Unknown option: %s', args{1});
    end    
  if isempty(args), break, end
end

% Normalize input to cell array gm_cell (nt x nK)
if isa(gm, 'gmdistribution')
  gm_cell = {gm};
elseif iscell(gm)
  gm_cell = gm;
else
  error('gm must be a gmdistribution or a cell array of gmdistribution objects.');
end

[nt, nK] = size(gm_cell);

% Initialize output
out = struct();
out.sort = doSort;
out.D = cell(nt, nK);
out.isort = cell(nt, nK);
out.minD = nan(nt, nK);
out.meanD = nan(nt, nK);

switch method
  case 'grid'
    [X,Y,Z] = ndgrid(vvec,vvec,vvec);
    XYZ = [X(:) Y(:) Z(:)];
end

%rmsF = cell(nt,nK);

for it = 1:nt
  for iK = 1:nK
    g = gm_cell{it,iK};    
    if isempty(g)
      continue
    end
    if ~isa(g, 'gmdistribution')
      error('gm{%d,%d} is not a gmdistribution.', it, iK);
    end
    K = g.NumComponents;
    
    % Sort by trace of Sigma (temperature , lowest first)
    % Not yet implemented
    if doSort
      isort = gmm_sort(g);
    else
      isort = 1:K;
    end
    out.isort{it,iK} = isort;

    rmsnorm = zeros(K,K);
    if strcmp(method,'montecarlo')
      xyz_tot = gmm_monte_carlo_sampling(g,N,1:K);
      [ptot, ~, ~] = gmm_evaluate_f(g,xyz_tot,1,'group',{1:K}); % evaluate f at xyz            
    end

    for ik1 = 1:K
      for ik2 = (ik1+1):K   
        kvec = [ik1 ik2];
        switch method
          case 'montecarlo' % generate points randomnly from the distribution
            g_merg = gmm_merge_components(g, {kvec});
            w_merg = sum(g.ComponentProportion(kvec));            
            xyz_p = gmm_monte_carlo_sampling(g,N,kvec); % generate points from the components specified in kvec
            %xyz_q = gmm_monte_carlo_sampling(g_merg.gmMerged,round(N*w_merg),1); % generate points from the components specified in kvec
            [~, ~, p] = gmm_evaluate_f(g,xyz_p,1,'group',{kvec}); % evaluate f at xyz
            [q, ~] = gmm_evaluate_f(g_merg.gmMerged,xyz_p,w_merg); % evaluate f at xyz
             %plot(p,q,'.'); grid on
            %%
            gmmFsummed = p;
            Fdiff = p-q;
          case 'grid'  % generate points on a grid
            g_merg = gmm_merge_components(g, {kvec});
            [g_Ftot, g_Fcomp, g_Fgrouped] = gmm_get_F(g,vvec,vvec,vvec,1,'group',{kvec}); 
            w = sum(g.ComponentProportion(kvec));
            [gmerg_Ftot, gmerg_Fcomp] = gmm_get_F(g_merg.gmMerged,vvec,vvec,vvec,w);
            
            % Caluclate difference (maxwellianity)
            gmmFsummed = squeeze(sum(g_Fcomp(:,:,:,kvec),4));
            Fdiff = gmerg_Ftot-gmmFsummed;
        end

        rmsnorm(ik1,ik2) = sqrt(sum(Fdiff.^2,'all'))/sqrt(sum(gmmFsummed.^2,'all'));
        %rmsnorm(ik1,ik2) = sum(abs(Fdiff),'all')/sum(abs(gmmFsummed),'all');
        rmsnorm(ik2,ik1) = rmsnorm(ik1,ik2);

        % do some awarding, populations with close drift speeds get lower
        % 'distance'.
        %dmu = norm(g_o.mu(ik1,:)*g_o.mu(ik2,:)');
        weight = sum(g.ComponentProportion([ik1, ik2]));
        %rmsnorm(ik1,ik2) = rmsnorm(ik1,ik2)*weight^0.5;
        %rmsnorm(ik2,ik1) = rmsnorm(ik2,ik1)*weight^0.5;

        
        if 0 % plot, debug, method=montecarlo
          %%
          vmax = max(sqrt(sum(XYZ(:).^2,2)));
          vmax = round(vmax*50)/50;
          vvec = -vmax:50:vmax;
          [pdf_tot,~,mid] = histcn(xyz_tot(:,3), vvec);
          [pdf_p,~,mid] = histcn(xyz_p(:,3), vvec);
          [pdf_q,~,mid] = histcn(xyz_q(:,3), vvec);
  
          figure; 
          hca = subplot(1,1,1);
                    
          plot(hca,mid{1},pdf_tot,mid{1},pdf_p,mid{1},pdf_q);
          legend(hca,{'Total, all components',sprintf('Summed components, ik = %g + %g',kvec(1),kvec(2)),sprintf('Merged, ik = %g + %g',kvec(1),kvec(2))})

          %legend(hca,{sprintf('ik = %g + %g',kvec(1),kvec(2)),})
          %legend(hca,{'Individual Gaussians summed','Merged into one Gaussian','Difference','Original components'},'box','off')
          %irf_legend(hca,{sprintf('ik = %g + %g',kvec(1),kvec(2)),sprintf('Sum(|Fmerge-Fsum|/|Fsum|)=%03.2f',rmsnorm(ik1,ik2))}',[0.02 0.98],'k')
        end
        if 0 % plot, debug, method=grid
          figure; 
          hca = subplot(1,1,1);
          plot(hca,vvec,squeeze(sum(gmmFsummed,[1 3])),vvec,squeeze(sum(gmerg_Ftot,[1 3])),vvec,squeeze(sum(Fdiff,[1 3])),vvec,squeeze(sum(g_Fcomp(:,:,:,kvec),[1 3])),'k--')
          %legend(hca,{sprintf('ik = %g + %g',kvec(1),kvec(2)),})
          legend(hca,{'Individual Gaussians summed','Merged into one Gaussian','Difference','Original components'},'box','off')
          irf_legend(hca,{sprintf('ik = %g + %g',kvec(1),kvec(2)),sprintf('Sum(|Fmerge-Fsum|/|Fsum|)=%03.2f',rmsnorm(ik1,ik2))}',[0.02 0.98],'k')
        end
      end
    end
    out.D{it,iK} = rmsnorm;

    tmp = rmsnorm;
    tmp(1:K+1:end) = NaN;
    out.minDB(it,iK) = min(tmp,[],'all','omitnan');
    out.meanDB(it,iK) = mean(tmp(~isnan(tmp)),'omitnan');
  end
end

% If input was a single gmdistribution, collapse outputs for convenience
if isa(gm, 'gmdistribution')
  out.D = out.D{1,1};
  out.isort = out.isort{1,1};
  out.minDB = out.minDB(1,1);
  out.meanDB = out.meanDB(1,1);
end
