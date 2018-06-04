%% Magnetic vector potential: Single time
% Az (Ay in GSM)
Az = @(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);
Bx = @(x,y,a,b,t) -b*y; % Bx = dyAz
By = @(x,y,a,b,t) -a*x; % By = -dxAz
Babs = @(x,y,a,b,t) sqrt(b^2*y.^2 + a^2*x.^2);

%divAz = @(x,y,a,b,t) (0+(a*x-b*y))-(t-1);
a = 1;
b = 7;
nA = 20;
t = 0;
nt = numel(t);



nx=251; ny=251;
xlim = 5;
ylim = 2;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);
[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));

% find areas of fluxtubes;

maxA = 14; % for a= 1, b = 7, t = 0;
minA = -14;
levA = minA:1:maxA;
nA = numel(levA);
A = Az(X(:,:),Y(:,:),a,b,t);
nMC = numel(A);
FVOL = nan(size(X));
fvol = nan(1,numel(levA));
for iA = 1:nA-1
  intA = [levA(iA) levA(iA+1)];
  indA = intersect(find(A>intA(1)),find(A<intA(2)));
  FVOL(indA) = numel(indA);
  fvol(iA) = numel(indA);
end

toplot=1:20:nx; 
BX = Bx(X,Y,a,b,t);
BY = By(X,Y,a,b,t);

% set up figure
fig = figure(10);
nrows=4; ncols=1; npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 1 % A
  hca = h(isub); isub = isub + 1;
  contourf(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),levA); 
  hold(hca,'on'); 
  quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end
if 1 % B abs
  hca = h(isub); isub = isub + 1;
  contourf(hca,X(:,:),Y(:,:),Babs(X(:,:),Y(:,:),a,b,t),minA:maxA); 
  hold(hca,'on'); 
  quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = '|B|';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end
if 1 % size of Fluxtube
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL);
  shading(hca,'flat');
  hold(hca,'on'); 
  contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  %quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  hca.CLim = [0 1];
  colormap('cool');
  colorbar('peer',hca)
end
if 1 % size of Fluxtube normalized
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL/max(max(FVOL)));
  shading(hca,'flat');
  hold(hca,'on'); 
  contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  %quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end

if 0 % divA
  hca = h(isub); isub = isub + 1;
  contourf(hca,X(:,:),Y(:,:),divAz(X(:,:),Y(:,:),a,b,t),nA); 
  %hold(hca,'on'); 
  %quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1];
  hold(hca,'off')
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end
%% Magnetic vector potential: Animation
xline=@(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;

nx=121; ny=121;
xlim = 4;
ylim = 2;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);

[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));
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
