ic = 1:4;
c_eval('mvaR?=irf.ts_vec_xyz(Rxyz?.time,[Rxyz?.data*v(1,:)'' Rxyz?.data*v(2,:)'' Rxyz?.data*v(3,:)'']);',ic)
tint_loc = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
c_eval('mvaR?loc = mvaR?.tlim(tint_loc);',ic)
c_eval('mvaR?loc = mean(mvaR?loc.data,1);',ic)
R0 = (mvaR1loc+mvaR2loc+mvaR3loc+mvaR4loc)/4;
c_eval('lmnR? = mvaR?loc-R0;',ic)

%
tint_plot = irf.tint('2015-10-16T10:33:23.00Z/2015-10-16T10:33:33.00Z');
Bl1 = [0; mvaB1.tlim(tint_plot).data(:,2);0;0];
Bl2 = [0; mvaB2.tlim(tint_plot).data(:,2);0;0];
Bl3 = [0; mvaB3.tlim(tint_plot).data(:,2);0;0];
Bl4 = [0; mvaB4.tlim(tint_plot).data(:,2);0;0];

ll1 = [0; linspace(0,100,numel(Bl1)-3)'; 100; 0];
ll2 = [0; linspace(0,100,numel(Bl2)-3)'; 100; 0];
ll3 = [0; linspace(0,100,numel(Bl3)-3)'; 100; 0];
ll4 = [0; linspace(0,100,numel(Bl4)-3)'; 100; 0];

patch(ll1+lmnR1(1),Bl1+lmnR1(2),Bl1*0+lmnR1(3),'black'); hold on;
patch(ll2+lmnR2(1),Bl2+lmnR2(2),Bl2*0+lmnR2(3),'red');
patch(ll3+lmnR3(1),Bl3+lmnR3(2),Bl3*0+lmnR3(3),'green');
patch(ll4+lmnR4(1),Bl4+lmnR4(2),Bl4*0+lmnR4(3),'blue');
hold off;
hca=gca;
hca.XLabel.String = 'L';
hca.YLabel.String = 'M';
hca.ZLabel.String = 'N';
