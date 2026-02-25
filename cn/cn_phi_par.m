function [phi] = cn_phi_par(E,B0,v)
% Find real E assuming it is parallel to B
B0=irf_resamp(B0,E);
B=c_coord_trans('gsm','dsi',B0,'cl_id',3);
%Bav=sum(B(:,2:4),1)/size(B,1);

theta=acos(B(3)./(B(:,2).^2+B(:,3).^2+B(:,4).^2));
Ez=sqrt(((E(:,2).^2+E(:,3).^2).*cos(theta).^2)./(sin(theta).^2));
Epar=[E(:,1) sqrt(E(:,2).^2+E(:,3).^2+Ez.^2)];


phi=irf_integrate(Epar)*v;
%irf_plot(phi)