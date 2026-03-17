function out = gmm_sort(gm,varargin)

var = squeeze(gm.Sigma(1,1,:) + gm.Sigma(2,2,:) + gm.Sigma(3,3,:));
[~,isort] = sort(var); % Ts
out = isort;