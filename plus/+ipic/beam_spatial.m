% Spatially localized beam that can be used to satisfy period boundary
% conditions for iPIC3D

%% a double gaussian (or some other greater power)
f = @(y,y0,J,l) J*exp(-(y-y0).^6/2/l/l);
L = 2; % width of the region contain ing one beam
y = linspace(-L,L,2000);
J = 1;
l = L/10;
plot(y,f(y,L/2,J,l),...
     y,f(y,-L/2,-J,l))
     %y,f(y,-L/2-l,-J,l)+f(y,-L/2+l,-J,l),...
     %y,f(y,L/2-1*l,J,l)+1*f(y,L/2+l,J,l)+f(y,L/2+3*l,J,l)+f(y,L/2-3*l,J,l)) 
     
%% Harris current sheet
L = 2; % width of the region contain ing one beam
y = linspace(-L,L,2000);
J = 1;
l = L/15;
B = @(y,l,B0) B0*tanh(y/l);
n = @(y,l,n0) n0*sech(y/l).^2;
n0 = 1; B0 = 1; y0 = L/2;
lw = 4*l;
plot(y,n(y,l,n0),y,B(y,l,n0),y,B(y,l,n0)-B(y-6*l,l,n0))
%plot(y,B(y-y0+lw,l,n0)-B(y-y0-lw,l,n0),...
%     y,-B(y+y0+lw,l,n0)+B(y+y0-lw,l,n0))
 
 %% With inspiration from Harris current sheet
L = 2; % width of the region contain ing one beam
y = linspace(-L,L,2000);
J = 1;
l = L/30;
J = @(y,l,J0) J0*tanh(y/l);
%J = @(y,l,lw,J0) J0*tanh((y-y0-lw)/l)-J0*tanh((y)-L);
B = @(y,l,J0) -(1-tanh(y/l).^2);
B0 = 1; 
y0 = L/2;
lw = 7*l;
Jtot = @(y,y0,l,lw,J0) J(y-y0+lw,l,J0)-J(y-y0-lw,l,J0)-J(y+y0+lw,l,J0)+J(y+y0-lw,l,J0);
Btot = @(y,y0,l,lw,B0) -B(y+y0+lw,l,J0)+B(y+y0-lw,l,J0)+B(y-y0+lw,l,J0)-B(y-y0-lw,l,J0);
if 0
plot(y,J(y-y0+lw,l,J0)-J(y-y0-lw,l,J0),...
     y,-J(y+y0+lw,l,J0)+J(y+y0-lw,l,J0),...
     y,-B(y+y0+lw,l,J0)+B(y+y0-lw,l,J0),...
     y,B(y-y0+lw,l,J0)-B(y-y0-lw,l,J0))
end
if 1
    %%
plot3(y,y*0,J(y-y0+lw,l,J0)-J(y-y0-lw,l,J0)-J(y+y0+lw,l,J0)+J(y+y0-lw,l,J0),...%y,y*0,-J(y+y0+lw,l,J0)+J(y+y0-lw,l,J0),...
    y,-B(y+y0+lw,l,J0)+B(y+y0-lw,l,J0)+B(y-y0+lw,l,J0)-B(y-y0-lw,l,J0),y*0)%,...
    %y,B(y-y0+lw,l,J0)-B(y-y0-lw,l,J0),y*0)
legend('J','B') 
xlabel('y')
ylabel('z')
zlabel('x')
end
if 1
    %%
    x = linspace(0,2*L,10);
    [X,Y] = meshgrid(x,y);    
    yind = 1:20:numel(y);
    pcolor(X,Y,Jtot(Y,y0,l,lw,J0)); shading flat; hold on;
    quiver(X(yind,:),Y(yind,:),-Btot(Y(yind,:),y0,l,lw,B0),Y(yind,:)*0,'k')
    %plot(y,Jtot(y,y0,l,lw,J0),y,Btot(y,y0,l,lw,B0))
    ylabel('z')
    xlabel('y')
    ch=colorbar;
    ylabel(ch,'J')
    legend(' ','B')
    colormap(cn.cmap('islands'))
    set(gca,'clim',2*[-1 1])
    title('Current and magnetic field')
    hold off
end
