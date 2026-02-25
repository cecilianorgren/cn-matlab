ux = 1.5;
uy = 2.1;
u = [ux uy];

vx = 1.2;
vy = 0.9;
v = [vx vy];


umag = sqrt(ux.^2 + uy.^2); 
uhat_x = ux/umag
uhat_y = uy/umag;
uhat = [uhat_x, uhat_y];

v_dot_uhat = vx*uhat_x + vy*uhat_y % v dot uhat

v_dot_uhat_x = uhat_x*v_dot_uhat
v_dot_uhat_y = uhat_y*v_dot_uhat;

v_minus_u_x = v_dot_uhat_x - ux;
v_minus_u_y = v_dot_uhat_y - uy;


hca = subplot(1,1,1);

%hq = gobjects(0);
quiver(hca,0,0,ux,uy,0);
hold(hca,'on')
quiver(hca,0,0,vx,vy,0)
quiver(hca,0,0,uhat_x,uhat_y,0)
quiver(hca,0,0,v_dot_uhat_x,v_dot_uhat_y,0)
quiver(hca,0,0,v_minus_u_x,v_minus_u_y,0)
hold(hca,'off')

hca.XLim = max(abs([hca.XLim hca.YLim]))*[-1 1];
hca.YLim = max(abs([hca.XLim hca.YLim]))*[-1 1];

axis square
%axis equal
grid on

legend({'u','v','uhat','v dot uhat','v - u'})

hq = findobj(gcf,'type','quiver'); hq = hq(end:-1:1);
for ii = 1:numel(hq)
  hq(ii).LineWidth = 2;
end

%%
N = 100;
xmin = -10; xmax = 10; x = linspace(xmin,xmax,N);
ymin = -10; ymax = 10; y = linspace(ymin,ymax,N);

ly = 1;
lx = 1;
Umax = 1;
fun_U = @(y,x) Umax*y.*exp(-(x/lx).^2-(y/ly).^2);

U = fun_U(X,Y);
contLevU = linspace(-Umax,Umax,50); 

nrows = 2;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,U)
shading(hca,'flat')
hold(hca,'on')
contour(hca,X,Y,U,contLevU,'k')
hold(hca,'off')
