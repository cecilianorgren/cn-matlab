

for kk=1:4; h(kk)=subplot(2,2,kk); end
isub=1;

hca=h(isub); isub=isub+1;
pcolor(hca,dg.R,dg.S,dg.wimax'/dg.wpi)
set(hca,'clim',[0 10]) 
colorbar('peer',hca)
hca=h(isub); isub=isub+1;
pcolor(hca,dg.R,dg.S,dg.wrmax'/dg.wpi)
set(hca,'clim',[0 50])
colorbar('peer',hca)

hca=h(isub); isub=isub+1;
pcolor(hca,cn.R,cn.S,cn.wimax/cn.wpi)
set(hca,'clim',[0 10]) 
colorbar('peer',hca)
hca=h(isub); isub=isub+1;
pcolor(hca,cn.R,cn.S,cn.wrmax/cn.wpi)
set(hca,'clim',[0 1000])
colorbar('peer',hca)
