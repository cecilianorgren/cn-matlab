syms vx vy vz Tx Ty Tz kB m vdx vdy vdz n pi

vtx = sqrt(2*kB*Tx/m);
vty = sqrt(2*kB*Ty/m);
vtz = sqrt(2*kB*Tz/m);

fx = n/((pi)^(1/2)*vtx)*exp(-(vx-vdx)^2/vtx/vtx);
fy = n/((pi)^(1/2)*vty)*exp(-(vy-vdy)^2/vty/vty);
fz = n/((pi)^(1/2)*vtz)*exp(-(vz-vdz)^2/vtz/vtz);

fxy = fx*fy;
fxyz = fx*fy*fz;

Tx = Ty;
vx = vy;
vdx = vdy;

simplify(fxy)

%% 
clear all
syms vx vy vz Tx Ty Tz kB m vdx vdy vdz n pi

vtx = sqrt(2*kB*Tx/m);
vty = sqrt(2*kB*Ty/m);
vtz = sqrt(2*kB*Tz/m);

fx = n/((pi)^(1/2)*vtx)*exp(-(vx-vdx)^2/vtx/vtx);
fy = n/((pi)^(1/2)*vty)*exp(-(vy-vdy)^2/vty/vty);
fz = n/((pi)^(1/2)*vtz)*exp(-(vz-vdz)^2/vtz/vtz);

fxy = fx*fy;
fxyz = fx*fy*fz;

Tx = Ty;
vx = vy;
vdx = vdy;

simplify(fxy)

%% 
clear all
syms vx vy vz Tx Ty Tz kB m vdx vdy vdz n pi

vtx = sqrt(2*kB*Tx/m);
vty = sqrt(2*kB*Ty/m);
vtz = sqrt(2*kB*Tz/m);

fx = n/((pi)^(1/2)*vtx)*exp(-(vx-vdx)^2/vtx/vtx);
fy = n/((pi)^(1/2)*vty)*exp(-(vy-vdy)^2/vty/vty);
fz = n/((pi)^(1/2)*vtz)*exp(-(vz-vdz)^2/vtz/vtz);

fxy = fx*fy;
fxyz = fx*fy*fz;

Tx = Ty;
vx = vy;
vdx = vdy;

simplify(fxy)