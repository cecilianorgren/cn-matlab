%% Based on simualtion
% turb = PIC('/Volumes/Fountain/Data/PIC/turbulencerun/data_h5/fields.h5');
doA = 1;

zlim = [-10 10];
twpe = 4900;
pic = turb.twpelim(twpe).zlim(zlim);
Bx = pic.Bx;
Bz = pic.Bz;
A = vector_potential(pic.xi,pic.zi,Bx,Bz);
vex = pic.vex;
vez = pic.vez;

hca = subplot(1,1,1);
imagesc(hca,pic.xi,pic.zi,smooth2(-vex,5,5)')
hca.CLim = [-4 4];
colormap(hca,pic_colors('blue_red'));

if doA
  hold(hca,'on')
  istep = 3;
  contour(hca,pic.xi(1:istep:end),pic.zi(1:istep:end),A(1:istep:end,1:istep:end)',-30:0.2:0,'k')
  hold(hca,'off')
end

hca.YDir = 'normal';
hca.XDir = 'normal';
hca.YLim = [-6 6];
hca.XLim = [20 47];


hcb = colorbar('peer',hca);
hcb.YTick = [];
hcb.YLabel.String = '  tailward            v_{ex}            Earthward';
hcb.FontSize = 11;
hca.XLabel.String = '<- Earthward                                        tailward ->';
hca.XTick = [];
hca.YTick = [];

%% Based on constructed vector potential

% topology 
% see also thesis.animation_reconnection

oline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2+b*(y-y0).^2);
oline   = @(x,y,x0,y0,a,b) 0.5*(exp(-a*(x-x0).^2-b*(y-y0).^2));
xline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2-b*(y-y0).^2);

%topology = @(x,y,x1,y1,x2,y2,a,b,c,d) xline(x,y,x1,y1,a,b) + oline(x,y,x2,y2,c,d);
topology = @(x,y) xline(x,y,0,0,1,7) + oline(x,y,3,0,2,2);
topology = @(x,y) 1*xline(x,y,0,0,1,7) + -1*oline(x,y,3,0,0.5,0.5);
topology = @(x,y) 1*xline(x,y,-5,0,1,7) + 1*xline(x,y,5,0,0.5,7);
topology = @(x,y) -1*(1*oline(x,y,-7,0,0.3,0.1) ...
                  + 1*oline(x,y,2,0,0.5,0.3) ...
                  + 1*0.004*xline(x,y,8,0,0.6,10));
%topology = @(x,y) oline(x,y,3,0,0.5,0.5);

% set up grid
nx=201; ny=201;
xlim = 12;
ylim = 8;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);
[X,Y]=meshgrid(x,y);
toplot=1:5:nx; 

% for time dependence, dont do this just yet
nt = 100;
t = linspace(1,10,nt);

% set up figure
figure(10) %figure('name','x-line')
ncols = 1;
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
Y1q = +4-x*4/40;
Bx_line = interp2(X,Y,Bx,X1q,Y1q);
By_line = interp2(X,Y,By,X1q,Y1q);
Babs_line = interp2(X,Y,Babs,X1q,Y1q);
vex = (topology(X(:,:),Y(:,:))-10).*By;
vex = (By).*(1-exp(-abs(Y)/5));

if 1 % A 
    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    %pcolor(hca,X(:,:),Y(:,:),vex); shading(hca,'flat')
   % hca.CLim = 0.02*[-1 1];
   % colorbar('peer',hca)
         
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),-2:0.1:2,'k');
    hold(hca,'on'); 
    contourf(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),[-2 0 2],'linewidth',1.5);
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),0*[1 1],'k','linewidth',1.5);
    colors = [ 0.99000    0.4    0.4; 0.2    0.8470    0.7410; 1 0.9 0; 0 0 0].^0.2;
    colormap(hca,colors)
    %colorbar
    
    
    %contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),-0.1025*[1 1],'k','linewidth',1.5);
    
    %axis(hca,'equal');  xlabel(hca,'X'); ylabel(hca,'Y');
    hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
    hca.YLim = [-6 6];
    hold(hca,'off')
    hca.XLabel.String = '<- Earthward                              tailward ->';
    %hca.YTick = [];
    %hca.XTick = [];
    hca.YLabel.String = '';
    hca.Box = 'on';
    hca.FontSize = 12;
    hca.Children = hca.Children(3:-1:1);
    
end
if 1 % |B|
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
if 0 % Bx
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
if 0 % By
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
if 0 % By
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
if 0 % By
    hca = h(isub); isub = isub + 1;
    plot(hca,x,[Bx_line; By_line; Babs_line])   
    hca.XLabel.String = 'X';
    hca.XLabel.String = 'B';
    legend(hca,'B_x','B_y','|B|','location','eastoutside')
end  


%% Based on constructed vector potential

% topology 
% see also thesis.animation_reconnection

oline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2+b*(y-y0).^2);
oline   = @(x,y,x0,y0,a,b) 0.5*(exp(-a*(x-x0).^2-b*(y-y0).^2));
xline   = @(x,y,x0,y0,a,b) 0.5*(a*(x-x0).^2-b*(y-y0).^2);

%topology = @(x,y,x1,y1,x2,y2,a,b,c,d) xline(x,y,x1,y1,a,b) + oline(x,y,x2,y2,c,d);
topology = @(x,y) xline(x,y,0,0,1,7) + oline(x,y,3,0,2,2);
topology = @(x,y) 1*xline(x,y,0,0,1,7) + -1*oline(x,y,3,0,0.5,0.5);
topology = @(x,y) 1*xline(x,y,-5,0,1,7) + 1*xline(x,y,5,0,0.5,7);
topology = @(x,y) -1*(1*oline(x,y,-3,0,0.1,0.2) ...
                  + 1*oline(x,y,-2,0,0.3,0.2) ...
                  + 1*0.02*xline(x,y,6,0,0.7,10));
%topology = @(x,y) oline(x,y,3,0,0.5,0.5);

topology_arr{1} = @(x,y) -1*(1*oline(x,y,-1,0,0.1,0.2) ...
                           + 1*oline(x,y,0,0,0.3,0.2) ...
                      + 1*0.02*xline(x,y,6,0,0.7,10));

topology_arr{2} = @(x,y) -1*(2*oline(x,y,-5,0,0.1,0.2) ...
                           + 1*oline(x,y,-10,0,0.3,0.2) ...
                      + 1*0.02*xline(x,y,6,0,0.7,10));

% set up grid
nx=201; ny=201;
xlim = 15;
ylim = 8;
x=linspace(-xlim,5,nx);
y=linspace(-ylim,ylim,ny);
[X,Y]=meshgrid(x,y);
toplot=1:5:nx; 

% for time dependence, dont do this just yet
nt = 100;
t = linspace(1,10,nt);

% set up figure
figure(10) %figure('name','x-line')
ncols = 1;
nrows = 1;
npanels = ncols*nrows;
for ip = 1:npanels
    h(ip) = subplot(nrows,ncols,ip);
end
isub = 1; 


%limA = [min(min(xline(X(:,:),Y(:,:),a,b))) max(max(xline(X(:,:),Y(:,:),a,b)))];
%limA(2) = limA(2)*2;
%linesA = linspace(limA(1),limA(2),20);
cmap = pic_colors('blue_red');
limA = [0 1];
Alev = -10:0.5:10;


if 1 % By (more than one time step)
  for itop = 1:numel(topology_arr)
    topology = topology_arr{itop};

    [dXx,dYx]=gradient(topology(X,Y));
    Bx = dYx;
    By = -dXx;
    Babs = sqrt(dXx.^2+dYx.^2);
    
    % interpolate to line
    X1q = x;
    Y1q = +4-x*4/40;
    Bx_line = interp2(X,Y,Bx,X1q,Y1q);
    By_line = interp2(X,Y,By,X1q,Y1q);
    Babs_line = interp2(X,Y,Babs,X1q,Y1q);
    vex = (topology(X(:,:),Y(:,:))-10).*By;
    vex = (By).*(1-exp(-abs(Y)/5));

    hca = h(isub); isub = isub + 1;
    %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
    pcolor(hca,X(:,:),Y(:,:),-By(:,:)); 
    shading(hca,'flat');
    hold(hca,'on')
    contour(hca,X(:,:),Y(:,:),topology(X(:,:),Y(:,:)),Alev+0.093,'k');     
    %quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    title(hca,'B_z')
    axis(hca,'equal');  
    xlabel(hca,'x'); 
    ylabel(hca,'z');
    %hca.XLim = xlim*[-1 1]; 
    %hca.YLim = ylim*[-1 1];
    hold(hca,'off')   
    hca.CLim = 0.05*[-1 1];
    colorbar('peer',hca)
    colormap(hca,cmap)
    %hca.XDir = 'reverse';
  end
end
