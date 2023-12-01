x = linspace(-3,3,10000);
z = linspace(-3,3,10000);

vx_x = @(x) x.*exp(-x.^2);
vy_x = @(x) -exp(-x.^2);
vxvy_x = @(x) vx_x(x).*vy_x(x);
vxdxvy = vx_x(x).*gradient(vy_x(x),x);

vz_z = @(z) -z.*exp(-z.^2);
vy_z = @(z) -exp(-z.^2);
vzvy_z = @(z) vz_z(z).*vy_z(z);
vzdzvy = vz_z(z).*gradient(vy_z(z),z);

hca = subplot(2,2,1);
plot(hca,x,vx_x(x),x,vy_x(x),x,vxvy_x(x))
hca.XLabel.String = 'x';
legend(hca,{'v_x','v_y','v_xv_y'},'location','best','box','off')

hca = subplot(2,2,3);
plot(hca,x,gradient(vxvy_x(x),x),x,vxdxvy)
hca.XLabel.String = 'x';
legend(hca,{'\partial_x(v_xv_y)','v_x\partial_xv_y'},'location','best','box','off')

hca = subplot(2,2,2);
plot(hca,z,vz_z(z),z,vy_z(z),z,vzvy_z(z))
hca.XLabel.String = 'z';
legend(hca,{'v_z','v_y','v_zv_y'},'location','best','box','off')

hca = subplot(2,2,4);
plot(hca,z,gradient(vzvy_z(z),z),z,vzdzvy)
hca.XLabel.String = 'z';
legend(hca,{'\partial_z(v_zv_y)','v_z\partial_zv_y'},'location','best','box','off')

h = findobj(gcf,'type','axes');
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 16;',1:numel(h))
c_eval('h(?).XLim = h(?).Children(1).XData([1 end]);',1:numel(h))

hl = findobj(h,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))