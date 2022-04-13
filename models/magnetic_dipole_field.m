function [Bx,By,Bz,R,EL,AZ,POL] = magnetic_dipole_field(x,y,z,k0)
% [Bx,By,Bz] = magnetic_dipole_field(x,y,z,k0)
%
% Spherical:
% B = r^(-3)*(3*(k_o*r)r-k0)
% Br = -2*0*r^(-3)*sin(lambda)
% Blambda = (k0/r)*cos(theta)



if numel(x) == 1
  r = sqrt(x.^2 + y.^2 + z.^2);
  polar = atand(sqrt(x.^2+y.^2)./z); % polar angle
  elev = -(polar-90);
  phi = atan2d(x,y);
  azim = phi;

  Br = -2*k0./(r.^3).*sin(elev);
  Btheta = k0./(r.^3).*cos(elev);
  Bphi = 0;

  Rot = [sin(polar)*cos(azim), cos(polar)*cos(azim), -sin(azim);...
       sin(polar)*sin(azim), cos(polar)*sin(azim),  cos(azim);...
       cos(polar),           -sin(polar),            0        ];
  Bxyz = Rot*[Br;Btheta;Bphi];

  Bx = Bxyz(1);
  By = Bxyz(2);
  Bz = Bxyz(3);
  R = r;
  EL = elev;
  AZ = azim;
elseif numel(x) > 1
  %sizedata = size(x);
  for ix = 1:size(x,1)
    for iy = 1:size(x,2)
      for iz = 1:size(x,3)
        
        r = sqrt(x(ix,iy,iz).^2 + y(ix,iy,iz).^2 + z(ix,iy,iz).^2);
        %polar = atand(sqrt(x(ix,iy,iz).^2+y(ix,iy,iz).^2)./z(ix,iy,iz)); % polar angle
        %elev = -(polar-90);
        
        %elev = atan(sqrt(x(ix,iy,iz).^2+y(ix,iy,iz).^2)./z(ix,iy,iz)); % elevation angle
        elev = atan(z(ix,iy,iz)./sqrt(x(ix,iy,iz).^2+y(ix,iy,iz).^2)); % elevation angle
        polar = pi/2 - elev;
        
        phi = atan2(x(ix,iy,iz),y(ix,iy,iz));
        azim = phi;

        Br = -2*k0./(r.^3).*sin(elev);
        Btheta = k0./(r.^3).*cos(elev);
        Bphi = 0;
        
        Rot = [sin(polar)*cos(azim), cos(polar)*cos(azim), -sin(azim);...
               sin(polar)*sin(azim), cos(polar)*sin(azim),  cos(azim);...
               cos(polar),           -sin(polar),            0        ];
     
        Bxyz = Rot*[Br; Btheta; Bphi];
        Bx(ix,iy,iz) = Bxyz(1);
        By(ix,iy,iz) = Bxyz(2);
        Bz(ix,iy,iz) = Bxyz(3);
        R(ix,iy,iz) = r;
        EL(ix,iy,iz) = elev;
        AZ(ix,iy,iz) = azim;
        POL(ix,iy,iz) = polar;
      end
    end
  end
end

