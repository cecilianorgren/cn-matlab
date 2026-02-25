function [y]=cn_xyz(inp,x,y,z)
% function output = cn_xyz(input,z,y,flag)
%
% input - input vector consists of columns [time in_x in_y in_z]
% x - 2nd direction
% z - for example magnetic field, [time bx by bz]
% output - output vector in new coordinates [time out_x out_y out_z] 

% [out]=irf_newxyz(inp,x,y,z)
% inp,out - column vector if more than 3 columns assume that first is time and 2nd-4th is X Y Z
% x,y,z - vectors (x=[xx xy xz], y= [yx yy yz]; )
%          that give new coordinates, if some is 0 then calculate it from other
% if x=0 -> x=yxz,z=xxy
% if y=0 -> y=zxx,x=yxz
% if z=0 -> z=xxy,y=zxx
%
% $Id$

if x==0
 x=irf_cross(y,z);
 z=irf_cross(x,y);
end
if y==0
 y=irf_cross(z,x);
 x=irf_cross(y,z); % to make sure x and z are orthogonal
end
if z==0
 z=irf_cross(x,y);
 y=irf_cross(z,x);
end

x=irf_resamp(irf_norm(x),inp);
y=irf_resamp(irf_norm(y),inp);
z=irf_resamp(irf_norm(z),inp);

out=inp;

if 0; size(out,2)==3,
  out(:,1)=inp(:,1:3)*x';
  out(:,2)=inp(:,1:3)*y';
  out(:,3)=inp(:,1:3)*z';
elseif 1; size(out,2)>3,
  out(:,2)=inp(:,2:4)*x(:,2:4)';
  out(:,3)=inp(:,2:4)*y(:,2:4)';
  out(:,4)=inp(:,2:4)*z(:,2:4)';
else
  irf_log('fcal','!!!!!! Coordinate transformation not possible !!!!!!!!!');
end


if 0
x=input;
if size(z,2)<4
    z=[input(1,1) z];
end
if size(z,2)<4
    x=[input(1,1) z];
end
  z=irf_resamp(z,x);

  nz=irf_norm(z); % along B
  ny=irf_norm(irf_cross(nz,irf_cross(y,nz))); % closest to given vn vector
  nx=irf_cross(ny,nz); % in (b x n) direction

% estimate e in new coordinates
  xx=irf_dot(x,nx,1);
  xy=irf_dot(x,ny,1);
  xz=irf_dot(x,nz,1);

  y=x; y(:,2)=xx;y(:,3)=xy;y(:,4)=xz;
end