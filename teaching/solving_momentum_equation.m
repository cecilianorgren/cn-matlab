syms vx(t) vy(t) vz(t) q m Bx By Bz Ex Ey Ez ax ay az t


ax = (q/m)*(Ex + vy*Bz-vz*By);
ay = (q/m)*(Ey + vz*Bx-vx*Bz);
az = (q/m)*(Ez + vx*By-vy*Bx);

eqn1 = diff(vx,t) == ax;
eqn2 = diff(vy,t) == ay;
eqn3 = diff(vz,t) == az;


[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [vx, vy, vz])
%%


X = linsolve(A,B);