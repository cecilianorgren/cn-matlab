b=[0 1 1];b=b/norm(b);
spaxis=[0 0 1];spaxis=spaxis/norm(spaxis);
sp=cn_cross(b,spaxis);sp=sp/norm(sp)
%vsp=cn_cross(dot(b,sp),b);
third=cn_cross(b,sp)

cn_plot3d([0,0,0],b,'b',sp,'sp',third,'third')
axis equal