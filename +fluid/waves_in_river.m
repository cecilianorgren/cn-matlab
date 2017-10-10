% Make explainatory figure for Ex 9.1



x = 0:3;
y = 0:0.1:200;

[X,Y] = meshgrid(x,y);

f_surf = @(x) cos(x);

Z = cos(Y);

surf(X,Y,cos(Y/10))
shading flat
axis off
set(gca,'clim',[-1.5 1.5])
view([3 1 2])
text(1.5, 100 ,1.6,'\leftarrow crest ','fontsize',30)
text(3, 58 ,-1.05,'trough \rightarrow','fontsize',30,'horizontalalignment','left')

set(gcf,'color','white'); % white background for figures (default is grey)
colormap('bone')

%%
figure(61)

quiver3(3,0,1,-1,0,0,1,'k')
hold on;
quiver3(3,0,1,0,0,1,1,'k')
quiver3(3,0,1,0,1,0,1,'k')
hold off;


%%

x = 0:2;
y = linspace(0,3,1000);

[X,Y] = meshgrid(x,y);

f_surf = @(x) cos(x);

Z = cos(Y);

surf(X,Y,cos(Y/0.13))
shading flat
axis off
set(gca,'clim',[-1.5 1.5])
view([3 1 2])
text(1.0, 1.35 ,1.8,'\leftarrow crest ','fontsize',30)
text(2, 0.6 ,-1.05,'trough \rightarrow','fontsize',30,'horizontalalignment','left')

set(gcf,'color','white'); % white background for figures (default is grey)
colormap('bone')
% add river bank
pX = [0 0 0 0];
pY = [0 3 3 0];
pZ = [-1 -1 1 1]*1.3;
hp=patch(pX,pY,pZ,'black');
bankcolor = [0.9064    0.918889    0.8889]*0.95;
hp.FaceColor = bankcolor;
hp.EdgeColor = bankcolor;
text(1.0, 2.5 ,2.7,'riverbank','fontsize',30)



%axis equal
hold on;
quiver3(2,0,1,-1,0,0,0.7,'k'); text(2, 0 ,2.3,'z ','fontsize',20,'verticalalignment','bottom')
quiver3(2,0,1,0,0,1,1.3,'k'); text(1.2, 0 ,1,'y ','fontsize',20,'verticalalignment','bottom')
quiver3(2,0,1,0,1,0,0.4,'k'); text(2, 0.4 ,1,'x ','fontsize',20,'verticalalignment','bottom')

% river flow direction
%quiver3(2,3,1,0,1,0,'k'); text(2, 0.4 ,1,'x ','fontsize',20,'verticalalignment','bottom')
%arrow([2 3 1],[0 1 0])
hold off;




