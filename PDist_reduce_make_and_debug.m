% PDist_reduce_make_and_debug

% Make reduced distribution
tintZoom = irf.tint('2017-07-06T13:53:50.00Z',25);
eint = [000 40000];
vint = [-Inf Inf];

iDist = iPDist1.tlim(tintZoom).elim(eint);
vi = gseVi1.tlim(iDist.time).resample(iDist);
iLine = dmpaB1.resample(iDist).norm;
iPlane1 = iLine.cross(irf.ts_vec_xyz(iLine.time,repmat([1 0 0],iLine.length,1)));
iPlane2 = iLine.cross(iPlane1);
if2D = iDist.reduce('2D',iPlane1,iPlane2,'vint',vint); % reduced distribution perp to B
%ef2D = eDist.reduce('2D',ePlane1,ePlane2,'vint',vint); % reduced distribution perp to B

%%
time = iDist(10).time;
h1 = subplot(2,1,1);
[a,b,h_all] = if2D.surf(h1,'tint',time,'printinfo');
axis(h1,'square')
shading(h1,'flat')

h2 = subplot(2,1,2);
xyz = [iPlane1.resample(time).data;iPlane2.resample(time).data;iLine.resample(time).data];
[hsf,plspec] = mms.plot_int_projection(h2,iDist,'t',time,'xyz',xyz,'vzint',vint,'colorbar',1);
%hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));

h1.XLim = h2.XLim;
h1.YLim = h2.YLim;
h1.CLim = h2.CLim;
h2.XGrid = 'on';
h2.YGrid = 'on';