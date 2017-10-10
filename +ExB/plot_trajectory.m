global lr lz phi0
close
modelnr=2;
ExB.model;
%%
% The electric field
Ex=@(x,y,z,lr,lz,phi0) (x/(lr*1e3)^2).*phi0.*exp(-0.5*(x/(lr*1e3)).^2-0.5*(y/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2);
Ey=@(x,y,z,lr,lz,phi0) (y/(lr*1e3)^2).*phi0.*exp(-0.5*(x/(lr*1e3)).^2-0.5*(y/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2);
Ez=@(x,y,z,lr,lz,phi0) (z/(lz*1e3)^2).*phi0.*exp(-0.5*(x/(lr*1e3)).^2-0.5*(y/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2);
%options = odeset('Events',@ExB.events,'InitialStep',2.5e-5);%,'OutputSel',1,'Refine',refine);

h = axes;
if 0 % Plot outer sphere
    [X,Y,Z] = ellipsoid(0,0,0,lr,lr,lz);
    C = 1e3*1e3*sqrt(Ex(X,Y,Z,lr,lz,phi0).^2+Ey(X,Y,Z,lr,lz,phi0).^2+Ez(X,Y,Z,lr,lz,phi0).^2);
    [posX] = find(X<-0.01);
    X(posX)=NaN; Y(posX)=NaN; Z(posX)=NaN; C(posX)=NaN;
    C = ones(size(Z,1)-1,size(Z,2)-1)*max(max(C));
    
    

    hs = surf(h,X,Y,Z,C); hold on;
    alpha(hs,0.2)
    %colorbar
end
if 0 % Plot inner sphere
    [X,Y,Z] = ellipsoid(0,0,0,lr*0.5,lr*0.5,lz*0.5);
    C = 1e3*1e3*sqrt(Ex(X,Y,Z,lr,lz,phi0).^2+Ey(X,Y,Z,lr,lz,phi0).^2+Ez(X,Y,Z,lr,lz,phi0).^2);
    [posX] = find(X<-0.01);
    X(posX)=NaN; Y(posX)=NaN; Z(posX)=NaN; C(posX)=NaN;
    C = ones(size(Z,1)-1,size(Z,2)-1)*max(max(C));
    hs = surf(h,X,Y,Z,C); hold on;
    alpha(hs,0.2)
    
end
if 1 % Plot inner sphere
    [X,Y,Z] = ellipsoid(0,0,0,lr,lr,lz);
    C = 1e3*1e3*sqrt(Ex(X,Y,Z,lr,lz,phi0).^2+Ey(X,Y,Z,lr,lz,phi0).^2);
    [posX] = find(X<-0.01);
    X(posX)=NaN; Y(posX)=NaN; Z(posX)=NaN; C(posX)=NaN;
    C = ones(size(Z,1)-1,size(Z,2)-1)*max(max(C));
    hs = surf(h,X,Y,Z,C); hold on;
    alpha(hs,0.2)
    
end
if 0 % Plot outer sphere
    [X,Y,Z] = ellipsoid(0,0,0,2*lr,2*lr,2*lz);
    C = 1e3*1e3*sqrt(Ex(X,Y,Z,lr,lz,phi0).^2+Ey(X,Y,Z,lr,lz,phi0).^2);
    [posX] = find(X<-0.01);
    X(posX)=NaN; Y(posX)=NaN; Z(posX)=NaN; C(posX)=NaN;
    C = ones(size(Z,1)-1,size(Z,2)-1)*max(max(C));
    hs = surf(h,X,Y,Z,C); hold on;
    alpha(hs,0.2)
    
end
if 0 % Plot outer sphere
    [X,Y,Z] = ellipsoid(0,0,0,lr*2,lr*2,lz*2);
    C = 1e3*1e3*sqrt(Ex(X,Y,Z,lr,lz,phi0).^2+Ey(X,Y,Z,lr,lz,phi0).^2+Ez(X,Y,Z,lr,lz,phi0).^2);
    [posX] = find(X<-0.01);
    X(posX)=NaN; Y(posX)=NaN; Z(posX)=NaN; C(posX)=NaN;
    C = ones(size(Z,1)-1,size(Z,2)-1)*max(max(C));
    hs = surf(h,X,Y,Z,C); hold on;
    alpha(hs,0.2)
    
end

%shading(h,'flat')
options = odeset('Events',@ExB.events,'InitialStep',2.5e-5);%,'OutputSel',1,'Refine',refine);
if 0 % Add a trapped particle
    tint = [0 0.05];
    x0 = 0;
    y0 = 0.8*lr*1e3;
    z0 = -lz*1e3;
    vtper = cn_eV2v(1600,'ev'); % eV -> km/s (perpendicular thermal velocity)
    vx0 = 3*vtper*1e3;
    vy0 = 0;
    vz0 = 1000*1e3;
    
    [t,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0],options);
    x = x_sol(2:end,1); 
    y = x_sol(2:end,2);
    z = x_sol(2:end,3);
    vx = x_sol(2:end,4);
    vy = x_sol(2:end,5);
    vz = x_sol(2:end,6);
    plot3(x*1e-3,y*1e-3,z*1e-3);
end

if 1 % Add a trapped particle
    tint = [0 0.05];
    x0 = 0;
    y0 = 0.5*lr*1e3;
    z0 = -0.2*lz*1e3;
    vtper = cn_eV2v(1600,'ev'); % eV -> km/s (perpendicular thermal velocity)
    vx0 = 0.5*vtper*1e3;
    vy0 = 0;
    vz0 = 2000*1e3;
    
    [t1,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0],options);
    x = x_sol(2:end,1); 
    y = x_sol(2:end,2);
    z = x_sol(2:end,3);
    vx = x_sol(2:end,4);
    vy = x_sol(2:end,5);
    vz = x_sol(2:end,6);
    plot3(x*1e-3,y*1e-3,z*1e-3,'r');
end

if 1 % Add a free particle
    tint = [0 0.05];
    x0 = 0;
    y0 = 1*lr*1e3;
    z0 = -3*lz*1e3;
    vtper = cn_eV2v(1600,'ev'); % eV -> km/s (perpendicular thermal velocity)
    vx0 = 1*vtper*1e3;
    vy0 = 0;
    vz0 = 2000*1e3;
    
    [t2,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0],options);
    x = x_sol(2:end,1); 
    y = x_sol(2:end,2);
    z = x_sol(2:end,3);
    vx = x_sol(2:end,4);
    vy = x_sol(2:end,5);
    vz = x_sol(2:end,6);
    plot3(x*1e-3,y*1e-3,z*1e-3,'r');
end

if 0 % Add a trapped particle
    tint = [0 0.05];
    x0 = 0;
    y0 = 0.2*lr*1e3;
    z0 = -0.2*lz*1e3;
    vtper = cn_eV2v(1600,'ev'); % eV -> km/s (perpendicular thermal velocity)
    vx0 = 2*vtper*1e3;
    vy0 = 0;
    vz0 = 2000*1e3;
    
    [t,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0],options);
    x = x_sol(2:end,1); 
    y = x_sol(2:end,2);
    z = x_sol(2:end,3);
    vx = x_sol(2:end,4);
    vy = x_sol(2:end,5);
    vz = x_sol(2:end,6);
    plot3(x*1e-3,y*1e-3,z*1e-3,'r');
end

if 1 % Add a trapped particle
    tint = [0 0.02];
    x0 = 2*lr*1e3;
    y0 = 0;
    z0 = -3*lz*1e3;
    vtper = cn_eV2v(1600,'ev'); % eV -> km/s (perpendicular thermal velocity)
    vx0 = 3*vtper*1e3;
    vy0 = 0;
    vz0 = 20000*1e3;
    
    [t3,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0],options);
    x = x_sol(2:end,1); 
    y = x_sol(2:end,2);
    z = x_sol(2:end,3);
    vx = x_sol(2:end,4);
    vy = x_sol(2:end,5);
    vz = x_sol(2:end,6);
    plot3(x*1e-3,y*1e-3,z*1e-3,'r');
end

grid off
axis equal
hold off;


