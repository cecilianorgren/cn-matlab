%% Rosetta
cd /Users/Cecilia/Research/Rosetta/

[V,F] = read_vertices_and_faces_from_obj_file('ESA_Rosetta_OSIRIS_67P_SHAP2P.obj');
V = V*2.5;
h = trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor',[0.7,0.7,0.7], 'EdgeColor', 'none');
light('Position',[-1.0,-1.0,100.0],'Style','infinite');
lighting phong;
axis equal;
grid off;
axis off;
set(gcf,'color','white');
view(-75,15)

colors = mms_colors('matlab');
h.FaceColor = colors(6,:).^0.8; mms_colors('4');

%% Cylinder
[X,Y,Z] = cylinder(zeros(2,1)+[0.1],10);
tri = delaunay(X(:),Y(:),Z(:));
%tr = triangulation(tri, X, Y, Z);
h = trisurf(tri,X,Y,Z,'FaceColor',[0.7,0.7,0.7], 'EdgeColor', 'none');
%h = trisurf(F,X,Y,Z,'FaceColor',[0.7,0.7,0.7], 'EdgeColor', 'none');
shading flat;
h.FaceAlpha = 0.5;
axis equal;
colors = mms_colors('matlab');
h.FaceColor = colors(6,:).^0.8; mms_colors('4');
light('Position',[-1.0,-1.0,100.0]);
lighting phong;
grid off;
axis off;
set(gcf,'color','white');
view(-75,15)

%% Cylinder
[X,Y,Z] = cylinder(zeros(10,1)+[0.1],100);
h = surf(X,Y,Z);
shading flat;
h.FaceAlpha = 0.5;
axis equal;
colors = mms_colors('matlab');
h.FaceColor = colors(6,:).^0.8.*[0.9 1 0.9]; mms_colors('4');
%light('Position',[-1.0,-1.0,100.0]);
grid off;
axis off;
set(gcf,'color','white');
view(-75,15)


%% 1/r^2 density gradient
fig = figure(17);
set(fig,'Position',[526   542   300   156])
r = 10:90;
n = 1./r.^2;
fontsize = 20;

h = subplot(1,1,1);
plot(h,r,n,'linewidth',2);
h.YLabel.String = 'n_n';
h.YTick = [];
h.XLabel.String = 'r';
h.XTick = [];
h.Box = 'off';
h.FontSize = fontsize;
h.XLim = [0 r(end)];

colors = mms_colors('matlab');
ht = text(20,0.5*h.YLim(2),'1/r^2','FontSize',fontsize,'Color',colors(1,:));

%%
% determine position of the axes
axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;
drawnow
% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);