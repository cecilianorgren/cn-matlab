function out = eh_costfunction_finite_tilt(params,x,y,z,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez)
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
lx = params(2);
ly = params(3);
lz = params(4);
angle = params(5);
x0 = params(6);
y0 = params(7);
z0 = params(8);
polar = params(9);

% Evaluate value at xyz-location
% For each search step, this evaluation is done for different phi0,lx0,
% etc.
for ic = 1:4
  ex(:,ic) = mf_Ex(x(:,ic),y(:,ic),z(:,ic),phi0,lx,ly,lz,angle,x0,y0,z0,polar);
  ey(:,ic) = mf_Ey(x(:,ic),y(:,ic),z(:,ic),phi0,lx,ly,lz,angle,x0,y0,z0,polar);
  ez(:,ic) = mf_Ez(x(:,ic),y(:,ic),z(:,ic),phi0,lx,ly,lz,angle,x0,y0,z0,polar);
end

% This should work for arrays to
weight_x = 1;
weight_y = 1;
weight_z = 1.0;
colors = mms_colors('1234');

if 0 % plot
  % Ex
  hca = subplot(3,1,1);
  plot(hca,z(:,1),ex(:,1),...
           z(:,2),ex(:,2),...
           z(:,3),ex(:,3),...
           z(:,4),ex(:,4))
  hold(hca,'on')
  plot(hca,z(:,1),Ex_data(:,1),...
           z(:,2),Ex_data(:,2),...
           z(:,3),Ex_data(:,3),...
           z(:,4),Ex_data(:,4))
  hold(hca,'off')
  
  % Ey
  hca = subplot(3,1,2);
  plot(hca,z(:,1),ey(:,1),...
           z(:,2),ey(:,2),...
           z(:,3),ey(:,3),...
           z(:,4),ey(:,4))
  hold(hca,'on')
  plot(hca,z(:,1),Ey_data(:,1),...
           z(:,2),Ey_data(:,2),...
           z(:,3),Ey_data(:,3),...
           z(:,4),Ey_data(:,4))
  hold(hca,'off')
  
  % Ez
  hca = subplot(3,1,3);
  plot(hca,z(:,1),ez(:,1),...
           z(:,2),ez(:,2),...
           z(:,3),ez(:,3),...
           z(:,4),ez(:,4))
  hold(hca,'on')
  plot(hca,z(:,1),Ez_data(:,1),...
           z(:,2),Ez_data(:,2),...
           z(:,3),Ez_data(:,3),...
           z(:,4),Ez_data(:,4))
  hold(hca,'off')
  1;%pause
elseif 0
  hca = subplot(3,1,1);
  plot(hca,ex)
  hold(hca,'on')
  plot(hca,Ex_data,'--')
  hold(hca,'off')
  
  hca = subplot(3,1,2);
  plot(hca,ey)
  hold(hca,'on')
  plot(hca,Ey_data,'--')
  hold(hca,'off')
  
  hca = subplot(3,1,3);
  plot(hca,ez)
  hold(hca,'on')
  plot(hca,Ez_data,'--')
  hold(hca,'off')
  pause(0.1)
end

out = sum(weight_x*(ex(:)-Ex_data(:)).^2 + weight_y*(ey(:)-Ey_data(:)).^2 + weight_z*(ez(:)-Ez_data(:)).^2);
%disp(sprintf('%g',out))

% function sse = sseval(x,tdata,ydata)
% These 'x' are the parameters I need to fit.
%   B0 = x(1);
%   B1 = x(2);
%   dt = x(3);
%   sse = sum((ydata - (B0 + B1*tanh(tdata/dt))).^2);
% end