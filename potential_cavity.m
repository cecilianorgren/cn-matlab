
lx = 15;
ly = 1;
V = @(x,y) tanh(x/lx).*exp(-y.^2/ly^2);
V = @(x,y) exp(-x.^2/lx^2).*exp(-y.^2/ly^2);

xvec = ly*linspace(-40,40,100);
yvec = ly*linspace(-3,3,100);
[X,Y] = ndgrid(xvec, yvec);

h = setup_subplots(2,1);
isub = 1;

ypick = ly*[0 0.2 0.4 0.6 0.8 1];
[Xp,Yp] = ndgrid(xvec, ypick);

hca = h(isub); isub = isub + 1;
contourf(hca,X,Y,V(X,Y))
hold(hca,'on')
plot(hca,Xp,Yp,'linewidth',1)
hold(hca,'off')
colormap(hca,flipdim(colormap('gray'),1))

hca = h(isub); isub = isub + 1;
plot(hca,xvec,V(Xp,Yp),'linewidth',1)