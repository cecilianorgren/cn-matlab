% topology 
% see also thesis.animation_reconnection

oline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2+b*(y-y0).^2);
oline   = @(x,y,x0,y0,a,b) 0.5*(exp(-a*(x-x0).^2-b*(y-y0).^2));
xline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2-b*(y-y0).^2);

%topology = @(x,y,x1,y1,x2,y2,a,b,c,d) xline(x,y,x1,y1,a,b) + oline(x,y,x2,y2,c,d);
topology = @(x,y) xline(x,y,0,0,1,7) + oline(x,y,3,0,2,2);
topology = @(x,y) 1*xline(x,y,0,0,1,7) + -1*oline(x,y,3,0,0.5,0.5);
topology = @(x,y) 1*xline(x,y,-5,0,1,7) + 1*xline(x,y,5,0,0.5,7);
topology = @(x,y) -1*(2*oline(x,y,-4,0,0.1,0.1) + 1*oline(x,y,1,0,0.2,0.2) + 1*0.005*xline(x,y,7,0,1,7));
%topology = @(x,y) oline(x,y,3,0,0.5,0.5);

% set up grid
nx=121; ny=121;
xlim = 10;
ylim = 5;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);
[X,Y]=meshgrid(x,y);
toplot=1:5:nx; 

% for time dependence, dont do this just yet
nt = 100;
t = linspace(1,10,nt);

% set up figure
figure(10) %figure('name','x-line')
ncols = 2;
nrows = 2;
npanels = ncols*nrows;
for ip = 1:npanels
    h(ip) = subplot(nrows,ncols,ip);
end
isub = 1; 


%limA = [min(min(xline(X(:,:),Y(:,:),a,b))) max(max(xline(X(:,:),Y(:,:),a,b)))];
%limA(2) = limA(2)*2;
%linesA = linspace(limA(1),limA(2),20);
limA = [0 1];

[dXx,dYx]=gradient(topology(X,Y));
Bx = dYx;
By = -dXx;
Babs = sqrt(dXx.^2+dYx.^2);

% interpolate to line
X1q = x;
Y1q = -0-x*2/20;
Bx_line = interpn(X',Y',Bx,X1q,Y1q);
By_line = interpn(X',Y',By,X1q,Y1q);

if 0 % A 
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    contourf(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),10); 
    hold(hca,'on'); 
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    title(hca,'B=(dA_z/dy,-dA_x/dx)')
    axis(hca,'equal');  xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hold(hca,'off')
    hca.CLim = limA;
    colorbar('peer',hca)
end
if 0 % |B|
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    pcolor(hca,X(:,:),Y(:,:),Babs(:,:)); 
    shading(hca,'flat');
    hold(hca,'on')
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),10,'k');     
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    title(hca,'|B|')
    axis(hca,'equal');  xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hold(hca,'off')   
    hca.CLim = 0.05*[-1 1];
    colorbar('peer',hca)
end
if 1 % Bx
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    pcolor(hca,X(:,:),Y(:,:),Bx(:,:)); 
    shading(hca,'flat');
    hold(hca,'on')
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),10,'k');     
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    axis(hca,'equal'); xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hold(hca,'off')   
    hca.CLim = 0.05*[-1 1];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_x';
end  
if 1 % By
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    pcolor(hca,X(:,:),Y(:,:),By(:,:)); 
    shading(hca,'flat');
    hold(hca,'on')
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),10,'k');     
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');    
    axis(hca,'equal');  xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hold(hca,'off')   
    hca.CLim = 0.05*[-1 1];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
end  
if 1 % By
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    pcolor(hca,X(:,:),Y(:,:),By(:,:)); 
    shading(hca,'flat');
    hold(hca,'on')
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),10,'k');     
    quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    %title(hca,'|B|')
    axis(hca,'equal');  xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hscat = scatter(hca,X1q,Y1q,By_line*0+10,By_line);
    %hscat.MarkerEdgeColor = [0 0 0];
    MarkerFaceColor = By_line;
    plot(hca,X1q,Y1q,'w')    
    hold(hca,'off')   
    hca.CLim = 0.05*[-1 1];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
end  
if 1 % By
    hca = h(isub); isub = isub + 1;
    plot(hca,x,[Bx_line; By_line])   
    hca.XLabel.String = 'X';
    hca.XLabel.String = 'B';
    legend(hca,'B_x','B_y','location','eastoutside')
end  
%%
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
