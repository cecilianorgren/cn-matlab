figure(36); posFig = [1 1 1676 955]; set(gcf,'position',posFig);
h = axes; set(h,'position',[0.0 0.0 1.0 1.0]);
scale = 0.1;
nx = fix(scale*posFig(4));
ny = fix(scale*posFig(3));
matrix = rand(nx,ny);
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[X,Y] = cn.meshgrid(x,y);

pcolor(h,matrix+X)
shading flat

box(h,'off')

colormap('winter')

