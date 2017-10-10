close
lz=20*1e3;
lr=10*1e3;

ellipse(lr*1e-3,lr*1e-3,0,0,0,[1 0 0],100,'xy')
ellipse(lr*1e-3,lz*1e-3,0,0,0,[1 0 0],100,'xz')
ellipse(lz*1e-3,lr*1e-3,0,0,0,[1 0 0],100,'yz')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
