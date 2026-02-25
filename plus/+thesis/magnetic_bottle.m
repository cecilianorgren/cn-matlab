% MR-kurs
%oline=@(x,y,a,b)(0.5*x.*y+0.5*(a*x.^2+b*y.^2));
%xline=@(x,y,a,b,t) (0.5*x.*y+0.5*(a*x.^2-b*y.^2))+t;
xline=@(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;

magnetic_bottle = @() z-z.^3;

nx=100; ny=100;
xlim = 4;
ylim = 2;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);

[X,Y]=meshgrid(x,y);
Z = magnetic_bottle(X);

contour(X,Y,Z)

%%
toplot=1:10:nx; 

nt = 100;
t = linspace(1,10,nt);

figure(10) %figure('name','x-line')
isub=1; npx=1; npy=1;

hca=subplot(npx,npy,isub); isub=isub+1;
a = 1;
b = 7;
limA = [min(min(xline(X(:,:),Y(:,:),a,b,t(nt)))) max(max(xline(X(:,:),Y(:,:),a,b,t(nt))))];
limA(2) = limA(2)*2;
linesA = linspace(limA(1),limA(2),20);
for it = 1:nt  
  [dXx,dYx]=gradient(xline(X,Y,a,b,t(it)));
  %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
  contourf(hca,X(:,:),Y(:,:),xline(X(:,:),Y(:,:),a,b,t(it)),linesA); 
  hold(hca,'on'); 
  quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
  title('B=(dA_z/dy,-dA_x/dx)')
  axis(hca,'equal');  xlabel('X'); ylabel('Y');
  hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
  hold(hca,'off')
  hca.CLim = limA;
  pause(0.1)
  %colorbar('peer',hca)
end
%%

figure('name','Az')
setupfigure
pl=[1 1 0]; 
isub=1; npx=sum(pl); npy=2;
if pl(1) % oline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    surf(hca,X,Y,oline(X,Y,1,1));   
    axis square; axis equal; 
    title('A_z=0.5*(x^2+y^2)');
    xlabel(hca,'X'); ylabel(hca,'Y'); zlabel(hca,'A_z'); 
    shading(hca,'flat')
end
if pl(2) % xline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    surf(hca,X,Y,xline(X,Y,1,1)); 
    axis square; axis equal; 
    title('A_z=0.5*(x^2-y^2)')
    xlabel(hca,'X'); ylabel(hca,'Y'); zlabel(hca,'A_z'); 
    shading(hca,'flat')
end
if pl(3) % xline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    surf(hca,X,Y,xline(X,Y,1,1)); 
    axis square; axis equal; 
    title('A_z=0.5*(x^2-2y^2)')
    xlabel(hca,'X'); ylabel(hca,'Y'); zlabel(hca,'A_z'); 
    shading(hca,'flat')
end
if pl(1) % xline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    contour(hca,X(toplot,toplot),Y(toplot,toplot),oline(X(toplot,toplot),Y(toplot,toplot),1,1)); 
    hold(hca,'on'); 
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYo(toplot,toplot),-dXo(toplot,toplot));
    title('B=(dA_z/dy,-dA_x/dx)')
    axis(hca,'square');  xlabel('X'); ylabel('Y');
    hold(hca,'off')
end
if pl(2) % xline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1)); 
    hold(hca,'on'); 
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot));
    title('B=(dA_z/dy,-dA_x/dx)')
    axis(hca,'square');  xlabel('X'); ylabel('Y');
    hold(hca,'off')
end
if pl(3) % xline Az
    hca=subplot(npx,npy,isub); isub=isub+1;
    contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,2)); 
    hold(hca,'on'); 
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx2(toplot,toplot),-dXx2(toplot,toplot));
    title('B=(dA_z/dy,-dA_x/dx)')
    axis(hca,'square');  xlabel('X'); ylabel('Y');
    hold(hca,'off')
end
