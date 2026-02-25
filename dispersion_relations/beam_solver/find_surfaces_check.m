% first run beam_solver_run_spis.m
% or load:
% 
% 
% 
N = numel(wiall(1,1,:,:,:,:));
N=N/2;
[nk ns nr nte2 nti ntries] = size(wiall);
nrows = 2;
ncols = round(N/nrows+0.4);
for kk = 1:N; h(kk) = subplot(nrows,ncols,kk); end
isub = 1;
ite1=1;
for ia = 1:nr
    for ib = 1:nte2
        for ic = 1:nti
            for id = 1:ntries                
                hca = h(isub); isub=isub+1;
                pcolor(hca,k,S,wiall(:,:,ia,ib,ic,id)')
                xlabel(hca,'k\lambda_{De}')
                ylabel(hca,'S')
                title(hca,{['nTe1=' num2str(Te1(ite1)) ', Te2=' num2str(Te2(ib)) ', Ti=' num2str(Ti(ic))],...
                            ['R=' num2str(R(ia)) ', try ' num2str(id)]})
                cb = colorbar('peer',hca);
                set(hca,'clim',0.3*[-1 1])                
            end
        end
    end
end
colormap(cn.cmap('bluered3'))