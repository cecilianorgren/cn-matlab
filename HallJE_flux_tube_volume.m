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



%nx=851; ny=551;
nx=121; ny=121;
xlim = 10;
ylim = 4;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);
[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));

% find areas of fluxtubes;

maxA = 20; % for a= 1, b = 7, t = 0;
minA = -20;
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

Bref = Babs(-4,2,a,b,t);
Thetaref = 90;
BABS = Babs(X(:,:),Y(:,:),a,b,t);
THETA = asind(sqrt(BABS/Bref*sind(Thetaref)));


toplot_x=1:20:nx; 
toplot_y=1:20:ny; 
BX = Bx(X,Y,a,b,t);
BY = By(X,Y,a,b,t);

%% figure
fig = figure(10);
nrows = 6; 
ncols = 1; 
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 1 % A
  hca = h(isub); isub = isub + 1;
  contourf(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),levA); 
  hold(hca,'on'); 
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
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
  pcolor(hca,X(:,:),Y(:,:),Babs(X(:,:),Y(:,:),a,b,t)); 
  shading(hca,'flat')
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  hc.LineWidth = 1.2;
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = '|B|';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  hca.CLim = [0 hca.CLim(2)];
  hca.CLim = [0 17];
  %pause(0.1)
  colorbar('peer',hca)
end
if 1 % adiabatic invariant, THETA
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),real(THETA)); 
  shading(hca,'flat')
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:2:maxA,'k');   
  %hc.LineWidth = 1.0;
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = '|B|';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = [0 hca.CLim(2)];
  hca.CLim = [0 90];
  %pause(0.1)
  colorbar('peer',hca)
end
if 1 % size of Fluxtube
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL);
  shading(hca,'flat')
  hold(hca,'on');  
  contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  colorbar('peer',hca)
end
if 1 % size of Fluxtube normalized
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL/max(max(FVOL)));
  shading(hca,'flat');
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  hc.LineWidth = 1.2;
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  hca.CLim = [0 1];
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end
if 1 % inverse size of fluxtube normalized
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),1./FVOL);
  shading(hca,'flat');
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  hc.LineWidth = 1.2;
  quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = [0 1];
  %hca.CLim = limA;
  %pause(0.1)
  colorbar('peer',hca)
end

colormap(cn.cmap('white_blue'));

for ipanel = 1:npanels
  %h(ipanel).XLim = [-5 5];
  %h(ipanel).YLim = [-2 2];
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
%% figure: illustrate Egedal2015
fig = figure(10);
nrows = 2; 
ncols = 1; 
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 0 % A
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
if 0 % B abs
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),Babs(X(:,:),Y(:,:),a,b,t)); 
  shading(hca,'flat')
  hold(hca,'on'); 
  contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');
  quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = '|B|';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  hca.CLim = [0 hca.CLim(2)];
  %pause(0.1)
  colorbar('peer',hca)
end
if 0 % size of Fluxtube
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL);
  shading(hca,'flat')
  hold(hca,'on');  
  contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  quiver(hca,X(toplot,toplot),Y(toplot,toplot),BX(toplot,toplot),BY(toplot,toplot),'k');
  hold(hca,'off'); 
  hca.Title.String = 'A_z';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  colorbar('peer',hca)
end
if 1 % size of Fluxtube normalized
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),FVOL/max(max(FVOL)));
  shading(hca,'flat');
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  hc.LineWidth = 1.2;
  hq = quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = 'Fluxtube volume';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  hca.CLim = [0 1];
  %hca.CLim = limA;
  %pause(0.1)
  hcb_tmp1 = colorbar('peer',hca)  
end
if 1 % inverse size of fluxtube normalized
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X(:,:),Y(:,:),1./FVOL);
  shading(hca,'flat');
  hold(hca,'on'); 
  [~,hc] = contour(hca,X(:,:),Y(:,:),Az(X(:,:),Y(:,:),a,b,t),minA:maxA,'k');   
  hc.LineWidth = 1.2;
  hq = quiver(hca,X(toplot_y,toplot_x),Y(toplot_y,toplot_x),BX(toplot_y,toplot_x),BY(toplot_y,toplot_x),'k');
  hold(hca,'off'); 
  hca.Title.String = 'n_e ~ 1 / Fluxtube volume';
  hca.XLabel.String = 'X'; 
  hca.YLabel.String = 'Y';   
  axis(hca,'equal');  
  hca.XLim = xlim*[-1 1]; 
  hca.YLim = ylim*[-1 1]; 
  hold(hca,'off')
  %hca.CLim = [0 1];
  %hca.CLim = limA;
  %pause(0.1)
  hca.CLim = [0 0.0005];
  hcb_tmp2 = colorbar('peer',hca);
end
colormap(cn.cmap('blue_white'));
colormap(cn.cmap('white_blue'));
%colormap('parula')
h(1).CLim = [0.35 1];
h(2).CLim = [0 0.0002];
for ipanel = 1:npanels
  h(ipanel).XLim = [-5 5];
  h(ipanel).YLim = [-2 2];
end

if 1 % formatting
  %hcb_tmp1.YLabel.String = h(1).Title.String;
  hcb_tmp1.Ticks = hcb_tmp1.YLim([1 end]);
  hcb_tmp1.TickLabels = {'small','large'};
  
  %hcb_tmp2.YLabel.String = h(2).Title.String;
  hcb_tmp2.Ticks = hcb_tmp2.YLim([1 end]);
  hcb_tmp2.TickLabels = {'small','large'};
  
  for ipanel = 1:npanels
    h(ipanel).XTick = [];
    h(ipanel).YTick = [];
    h(ipanel).XLabel.String = [];
    h(ipanel).YLabel.String = [];    
    h(ipanel).Title.FontSize = 14;    
    
  end

end

%% Magnetic vector potential: Animation
xline=@(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline=@(x,y,a,b,t) ((t-1)/20).*0.6*a*x.^2 -((t-1)/10)*2*b*y.^2 - 1*(t-10);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline=@(x,y,a,b,t) t*b*y.^2;%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;

nx=121; ny=121;
xlim = 4;
ylim = 2;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);

[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));
toplot=1:10:nx; 

nt = 40;
t = linspace(1,20,nt);

figure(10) %figure('name','x-line')
isub=1; npx=1; npy=1;

hca=subplot(npx,npy,isub); isub=isub+1;
a = 1;
b = 7;
limA = [min(min(xline(X(:,:),Y(:,:),a,b,t(nt)))) max(max(xline(X(:,:),Y(:,:),a,b,t(nt))))];
limA(2) = limA(2)*2;
linesA = 1*linspace(limA(1),limA(2),20);
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

%% Magnetic vector potential: Animation
xline=@(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline=@(x,y,a,b,t) 0+abs(0.0*(a*y))+(0+0.5*(0*a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline=@(x,y,a,b,t) 0-abs(abs(a*y).^(0.01*t))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline=@(x,y,a,b,t) abs(a*y).^(1-0.8*t);

nx=121; ny=121;
xlim = 2;
ylim = 1;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);

[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));
toplot=1:10:nx; 

nt = 10;
t = linspace(0,1,nt);

figure(10) %figure('name','x-line')
isub=1; npx=1; npy=1;

hca=subplot(npx,npy,isub); isub=isub+1;
a = 1;
b = 7;
limA = [min(min(xline(X(:,:),Y(:,:),a,b,t(nt)))) max(max(xline(X(:,:),Y(:,:),a,b,t(nt))))];
%limA(2) = limA(2)*2;
linesA = linspace(limA(1),limA(2),20);
linesA = 20;
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
  %hca.CLim = limA;
  pause(0.1)
  %colorbar('peer',hca)
end

%% Magnetic vector potential: Animation
tend = 20;

xline=@(x,y,a,b,t) (0+0.5*(a*x.^2-b*y.^2))-(t-1);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
%xline=@(x,y,a,b,t) ((t-1)/20).*0.6*a*x.^2 -((t-1)/10)*2*b*y.^2 - 1*(t-10);%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
xline_pretend = @(x,y,a,b,t) -t*b*y.^2;%+t*0.05*(a*x.^2+b*y.^2);%-0.5*(a*x.^2+b*y.^2)*0.02*t ;%-t-(t-5)*1.1*b*y.^2;
%xline=@(x,y,a,b,t) -1*tend*b*y.^2 + (t/tend)*0.2*((0+0.5*(a*x.^2-0*b*(y-0*tend*y).^2))-10*(t-0));
xline_posttend = @(x,y,a,b,t) -1*tend*b*y.^2 + (t-tend)*0.5*a*x.^2-10*(t-tend);

nx=121; ny=121;
xlim = 4;
ylim = 2;
x=linspace(-xlim,xlim,nx);
y=linspace(-ylim,ylim,ny);

[X,Y]=meshgrid(x,y);
%[dXo,dYo]=gradient(oline(X,Y,1,1));
%[dXx2,dYx2]=gradient(xline(X,Y,1,2));
toplot=1:10:nx; 

doQuivers = 0;
doGif = 1;
nt = 40;

t = linspace(0,70,nt);

colors = mms_colors('matlab');
figure(10) %figure('name','x-line')
isub=1; npx=1; npy=1;

hca=subplot(npx,npy,isub); isub=isub+1;
hca.Position = [0 0 1 1];
a = 1;
b = 10;
%limA = [min(min(xline(X(:,:),Y(:,:),a,b,t(nt)))) max(max(xline(X(:,:),Y(:,:),a,b,t(nt))))];
%limA(2) = limA(2)*2;
  limA = [-2000 -100];
linesA = 1*linspace(limA(1),limA(2),20);
for it = 1:nt  
  if t(it)<tend
    xline = xline_pretend;
  else
    xline = xline_posttend;
  end
  [dXx,dYx]=gradient(xline(X,Y,a,b,t(it)));
  %contour(hca,X(toplot,toplot),Y(toplot,toplot),xline(X(toplot,toplot),Y(toplot,toplot),1,1,t(it))); 
  [cc,hc] = contourf(hca,X(:,:),Y(:,:),xline(X(:,:),Y(:,:),a,b,t(it)),linesA); 
  hc.LineWidth = 2;
  hc.Color = colors(1,:);
  if doQuivers
    hold(hca,'on'); 
    hq = quiver(hca,X(toplot,toplot),Y(toplot,toplot),dYx(toplot,toplot),-dXx(toplot,toplot),'k');
    hq.Color = colors(1,:);
    hq.LineWidth = 1.5;
    hold(hca,'off')
  end
  %title('B=(dA_z/dy,-dA_x/dx)')
  axis(hca,'equal');  xlabel('X'); ylabel('Y');
  hca.XLim = xlim*[-1 1]; hca.YLim = ylim*[-1 1];
  hca.CLim = limA;
  %colorbar
  %colormap(cn.cmap('blue_white'))
  colormap([0 0 0; 0 0 0])
  axis(hca,'off')
  pause(0.1)
  
  %colorbar('peer',hca)
  i_frame = it;
  if doGif
    drawnow;    
    f = getframe(gcf);
    %A(i_frame) = f;
    if it == 1 % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,20) = 0;
    else
        im(:,:,1,i_frame) = rgb2ind(f.cdata,map,'nodither');
    end
  end  
end

if 0*doGif
  %movie(fg,A,20);
  imwrite(im,map,sprintf('%sgifs/downsampling_t_step%g_reverse_1.gif',eventPath,t_step),'DelayTime',0.1,'LoopCount',0)
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
