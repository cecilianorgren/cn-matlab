tintJ = irf.tint('2015-11-12T07:19:15.00Z/2015-11-12T07:19:45.00Z');

J = gseJ1.tlim(tintJ);
cumJ = cumsum(J.data,1);
v = 50;
x = [v*(J.time-J.time.start)];
dx = (x(2)-x(1))*1e3;  % m
xyz = [x x*0 x*0];

hca = subplot(2,1,1);
plot_quivers(hca,J,xyz,'k');

hca = subplot(2,1,2);
plot_quivers(hca,cumJ(1:10:end,:)*dx,xyz(1:10:end,:),'k');
maxlim = max([abs(hca.YLim) abs(hca.ZLim)]);
hca.YLim = maxlim*[-1 1];
hca.ZLim = maxlim*[-1 1];