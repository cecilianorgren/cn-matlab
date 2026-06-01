function out = gmm_compare_vdfs(g1,g1)
% Idea is to compare the overall end gaussianity, or difference between
% original and merged gmm, to see how many groups are needed for a given
% end threshold.

[~, ~, p_] = gmm_evaluate_f(g_o,xyz_p_o,1,'group',{[1 2]});
[q_, ~] = gmm_evaluate_f(g_m,xyz_p_o,1); % evaluate f at xyz