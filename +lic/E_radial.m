E = @(r,lr) ((r/lr^2).*exp(-r.^2/2/(lr^2)));

nx=400;
ny=nx;
rmax=4;
lr=1;
x = linspace(-rmax,rmax,nx)*lr;
y = x;

set(gcf,'defaultTextFontSize',16);
set(gcf,'defaultAxesFontSize',16);


[X,Y]=meshgrid(x,y);

R = sqrt(X.^2+Y.^2);

h=axes;
pcolor(h,X,Y,E(R,lr)')
shading flat
axis square
hc=colorbar('peer',h);
%set(gcf,'defaultTextFontSize',16);
%set(gcf,'defaultAxesFontSize',16);
xlabel(h,'x/l_r')
ylabel(h,'y/l_r')
ylabel(hc,'E  [arb. units]')
colormap(cmap)

%%
if 0
%%
set(gcf,'defaultTextFontSize',16);
set(gcf,'defaultAxesFontSize',16);
plot(0,1)
ylabel('\rho _e/l_r','fontsize',20)
xlabel('r_{gc}/l_r','fontsize',20)
hc=colorbar;
ylabel(hc,'sign(v'')\times log_{10}|v''|','fontsize',20)
    
end