function gm_out = gmm_sort_components(gm)

sort_var = squeeze((gm.Sigma(1,1,:)+gm.Sigma(2,2,:)+gm.Sigma(3,3,:))/3);
[~,isort] = sort(sort_var); % Ts
        
mu = gm.mu(isort,:);
Sigma = gm.Sigma(:,:,isort);
ComponentProportion = gm.ComponentProportion(isort);
gm_out = gm;
%gm_out.mu = mu;
%gm_out.Sigma = Sigma;
%gm_out.ComponentProportion = ComponentProportion;

gm_out = gmdistribution(mu, Sigma, ComponentProportion);
        