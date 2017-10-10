% magnetic dipole field in cylindrical coordinates
B=@(z,rho)((z.^2+rho.^2).^(-3/2).*sqrt(1+3*z.^2./(z.^2+rho.^2)));

nPlot=4;
for ii=1:nPlot; h(ii)=subplot(nPlot,1,ii); end
iSub=1;

if 1
    hca=h(iSub); iSub=iSub+1;    
    z=linspace(-4,4,100); rho=0.1;
    plot(hca,z,B(z,rho)); 
    ylabel(hca,'B'); xlabel(hca,'z'); 
    title(hca,['\rho=',num2str(rho)]);
end
if 1
    hca=h(iSub); iSub=iSub+1;    
    z=linspace(-20,20,100); rho=5;
    plot(hca,z,B(z,rho)); 
    ylabel(hca,'B'); xlabel(hca,'z'); 
    title(hca,['\rho=',num2str(rho)]);
end
if 1
    hca=h(iSub); iSub=iSub+1;    
    z=0.1; rho=linspace(-1,1,100);
    plot(hca,rho,B(z,rho)); 
    ylabel(hca,'B'); xlabel(hca,'rho'); 
    title(hca,['z=',num2str(z)]);
end
if 1
    hca=h(iSub); iSub=iSub+1;    
    z=logspace(-3,0,400); rho=logspace(-1,0,400);
    [RHO,Z]=meshgrid(rho,z);
    surf(hca,RHO,Z,B(Z,RHO)); 
    ylabel(hca,'z'); xlabel(hca,'rho'); 
    set(hca,'zscale','log')
    shading(hca,'flat');
    hca_c=colorbar;
    caxis([0 200]);
    view(hca,[0 0 1])
end