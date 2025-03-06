function x_res  = eom_pic_external_griddedInterpolant(t,x_vect,fEx,fEy,fEz,fBx,fBy,fBz,m,q)
  x = x_vect(1);
  y = x_vect(2);
  z = x_vect(3);
  vx = x_vect(4);
  vy = x_vect(5);
  vz = x_vect(6);

  if isnan(x)
    1;
  end
  %disp(sprintf('t = %g, x = %g, y = %g, z = %g',t,x,y,z))

  %method = 'spline';        

  Ex = fEx(x,z);
  Ey = fEy(x,z);
  Ez = fEz(x,z);
  Bx = fBx(x,z);
  By = fBy(x,z);
  Bz = fBz(x,z);
  
  % Equations to be solved
  x_res = zeros(6,1);
  x_res(1) = vx; % dx/dt = vx;
  x_res(2) = vy; % dy/dt = vy;
  x_res(3) = vz; % dz/dt = vz;
  x_res(4) = (q/m)*(Ex + vy*Bz - vz*By);
  x_res(5) = (q/m)*(Ey + vz*Bx - vx*Bz);
  x_res(6) = (q/m)*(Ez + vx*By - vy*Bx);                                              

end 