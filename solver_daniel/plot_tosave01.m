index = 0;
eval(['toplot=tosave' num2str(index) ';'])
load
for kk=1:4; h(kk)=subplot(2,2,kk); end
isub=1;

hca=h(isub); isub=isub+1;
pcolor(hca,R,S,toplot.wrmax'/wpi)
set(hca,'clim',[0 20]) 
colorbar('peer',hca)
hca=h(isub); isub=isub+1;
pcolor(hca,R,S,toplot.wimax'/wpi)
set(hca,'clim',[0 5])
colorbar('peer',hca)
hca=h(isub); isub=isub+1;
pcolor(hca,R,S,toplot.vphmax'/veth1)
set(hca,'clim',[0 1.8])
colorbar('peer',hca)
hca=h(isub); isub=isub+1;
pcolor(hca,R,S,toplot.kmax'*Ld)
set(hca,'clim',[0 2])
colorbar('peer',hca)
if 0
hca=h(isub); isub=isub+1;
pcolor(hca,R,S,toplot.residual/wpi)
set(hca,'clim',[0 1.2])
colorbar('peer',hca)
end
colormap(cn.cmap('bluered3'))

%%
index = 0;
eval(['toplot=tosave' num2str(index) ';'])

    iR = 3;
    wr_kS = nan(nk,nv);
    wi_kS = nan(nk,nv);
    for iv = 1:nv;
        wr_kS(1:size(wr{iR,iv},2),iv) = toplot.wr{iR,iv};
        wi_kS(1:size(wr{iR,iv},2),iv) = toplot.wi{iR,iv};
        res_kS(1:size(wr{iR,iv},2),iv) = toplot.residual{iR,iv};
    end
    subplot(3,1,1)
    pcolor(kvec*Ld,S,wr_kS'/wpi)
    set(gca,'clim',[0 30])
    shading('flat')
    colorbar
    subplot(3,1,2)
    pcolor(kvec*Ld,S,wi_kS'/wpi)
    set(gca,'clim',[0 30])
    shading('flat')
    colorbar
    subplot(3,1,3)
    pcolor(kvec*Ld,S,log10(res_kS'/wpi))    
    %set(gca,'clim',[0 ]))
    shading('flat')
    cb = colorbar;
    set(cb,'ytick',[-10:10])
    ylabel(cb,'log10(residual)')
    
%%
    iR = 3;
    wr_kS = nan(nk,nv);
    wi_kS = nan(nk,nv);
    for iv = 1:nv;
        wr_kS(1:size(wr{iR,iv},2),iv) = wr{iR,iv};
        wi_kS(1:size(wr{iR,iv},2),iv) = wi{iR,iv};
        res_kS(1:size(wr{iR,iv},2),iv) = residual{iR,iv};
    end
    subplot(3,1,1)    
    pcolor(kvec*Ld,S,wr_kS'/wpi)
    title(['R=' num2str(R(iR))])
    set(gca,'clim',[0 30])
    shading('flat')
    colorbar
    subplot(3,1,2)
    pcolor(kvec*Ld,S,wi_kS'/wpi)
    set(gca,'clim',8*[-1 1])
    shading('flat')
    cb=colorbar;
    set(cb,'ylim',8*[-1 1])
    subplot(3,1,3)
    pcolor(kvec*Ld,S,log10(res_kS'/wpi))    
    %set(gca,'clim',[0 ]))
    shading('flat')
    cb = colorbar;
    set(cb,'ytick',[-10:10])
    ylabel(cb,'log10(residual)')