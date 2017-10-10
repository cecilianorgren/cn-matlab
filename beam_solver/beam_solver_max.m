function [k_max,wr_max,wi_max,vph_max,vph] = beam_solver_max(k,wr,wi,doPlot)
% take wi and wr from beam_solver and kind maximum growth rates, constructs
% phase velocities.
%
% Assume input omega is normalized to electron plasma frequency and that k
% is normalized to 1/l_De.

%doPlot = 0;

wr_orig = wr;
wi_orig = wi;
k_orig = k;

nk = size(wr,1);
nv = size(wr,2);
nR = size(wr,3); 
nTe2 = size(wr,4);
nTi = size(wr,5);

% make phase velocity
k_mat = repmat(k',1,nv,nR,nTe2,nTi);
vph = wr./k_mat/sqrt(2); % phase velocity in terms of background thermal speed

%% clear up the matrix from phoney values
% take away solutions that have negative frequencies
wi(wr_orig<0) = NaN;
vph(wr_orig<0) = NaN;
wr(wr_orig<0) = NaN;

if 0 % take away solutions that have too high frequencies
    wr_lim = 2.4;
    wi(wr_orig>wr_lim) = NaN;
    vph(wr_orig>wr_lim) = NaN;
    wr(wr_orig>wr_lim) = NaN;
end
%% find the point where the growth is max
% find the max for each R and S, i.e. one per k-vector
k_max = zeros(nv,nR,nTe2,nTi);
k_max_ind = zeros(nv,nR,nTe2,nTi);
wr_max = zeros(nv,nR,nTe2,nTi);
wi_max = zeros(nv,nR,nTe2,nTi);
vph_max = zeros(nv,nR,nTe2,nTi);

for iTi = 1:nTi;
    for iTe2 = 1:nTe2;        
        for iR = 1:nR
            for iS = 1:nv
                [ind,~] = find(wi(:,iS,iR,iTe2,iTi)==max(wi(:,iS,iR,iTe2,iTi)));                
                if ~isempty(ind)
                    if numel(ind)>1
                        disp('numerous max, using first')
                        ind=ind(1);
                    end
                    try
                        disp(['iR=' num2str(iR) ', iS=' num2str(iS) ', max_ind=' num2str(ind) ', found max'])
                    catch
                        1
                    end
                    
                    k_max_ind(iS,iR,iTe2,iTi) = ind;                    
                    vph_max(iS,iR,iTe2,iTi) = vph(ind,iS,iR,iTe2,iTi);
                    k_max(iS,iR,iTe2,iTi) = k_mat(ind,iS,iR,iTe2,iTi);
                    wi_max(iS,iR,iTe2,iTi) = wi(ind,iS,iR,iTe2,iTi);
                    wr_max(iS,iR,iTe2,iTi) = wr(ind,iS,iR,iTe2,iTi);
                else
                    disp(['iR=' num2str(iR) ', iS=' num2str(iS) ', max_ind=' num2str(ind) ', empty ind'])                  
                    vph_max(iS,iR,iTe2,iTi) = NaN;
                    k_mat(iS,iR,iTe2,iTi) = NaN;
                    wi_max(iS,iR,iTe2,iTi) = NaN;
                    wr_max(iS,iR,iTe2,iTi) = NaN;
                end
                if doPlot
                    plot(k_orig,wr(:,iS,iR,iTe2,iTi))
                    %pause
                end
                %disp(['iR=' num2str(iR) ', iS=' num2str(iS) ', max_ind=' num2str(ind) ', vph_max_ind=' num2str(vph_max_store(iS,iR))])        
            end
        end
    end
end