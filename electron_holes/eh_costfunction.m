function out = eh_costfunction(params,x,y,z,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez)
% x, y, z are predefined from tetrahedron configuration
% x0, y0, z0 displaces it if needed
% Cost function needs to be made for all four spacecraft. 
% Give Ex,y,z_data as a 4x1 array from observations, same with x,y,z: x,y,z
% is the relative location of the spacecraft.
% To get the final cost function, we need to supply the data from the
% observations such that only params is the input:
% 
% xdata = [R1(1),R2(1),R3(1),R4(1)];
% ydata = [R1(2),R2(2),R3(2),R4(2)];
% zdata = [R1(3),R2(3),R3(3),R4(3)];
% Ex_data = [E1(1),E2(1),E3(1),E4(1)];
% Ey_data = [E1(2),E2(2),E3(2),E4(2)];
% Ez_data = [E1(3),E2(3),E3(3),E4(3)];
% cost_function = @(params) eh_costfunction(params,xdata,ydata,zdata,mf_Ex,mf_Ey,mf_Ez)
%
% Alternatively, we can also pass xdata,ydata,zdata into mf_Ex,y,z...


phi0 = params(1);
lx0 = params(2);
ly0 = params(3);
lz0 = params(4);
angle = params(5);
x0 = params(6);
y0 = params(7);
z0 = params(8);

% Evaluate value at xyz-location
% For each search step, this evaluation is done for different phi0,lx0,
% etc.
for ic = 1:4
  ex(ic) = mf_Ex(x(ic),y(ic),z(ic),phi0,lx0,ly0,lz0,angle,x0,y0,z0);
  ey(ic) = mf_Ey(x(ic),y(ic),z(ic),phi0,lx0,ly0,lz0,angle,x0,y0,z0);
  ez(ic) = mf_Ez(x(ic),y(ic),z(ic),phi0,lx0,ly0,lz0,angle,x0,y0,z0);
end

% This should work for arrays to
out = sum((ex-Ex_data).^2 + (ey-Ey_data).^2 + (ez-Ez_data).^2);

if params(1) < 0 % will not allow negative density
  out = inf;
end
% function sse = sseval(x,tdata,ydata)
% These 'x' are the parameters I need to fit.
%   B0 = x(1);
%   B1 = x(2);
%   dt = x(3);
%   sse = sum((ydata - (B0 + B1*tanh(tdata/dt))).^2);
% end