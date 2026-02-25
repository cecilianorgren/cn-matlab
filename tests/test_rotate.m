phi=30;
theta=0;

o=[0 0 0];
v=[1 0 0];
rotate(v,[theta phi],1)
quiver3(o(1),o(2),o(3),v(1),v(2),v(3)); hold on;
quiver3(o(1),o(2),o(3),v(1),v(2),v(3))


