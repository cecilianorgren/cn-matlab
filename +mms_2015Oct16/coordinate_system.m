%% LMN

% N
tint = irf.tint('2015-10-16T10:33:22.000Z/2015-10-16T10:33:32.000Z');
[out,l,v] = irf_minvar(jbrst.tlim(tint));
lmN = v(3,:);

tint = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.400Z');
[out,l,v] = irf_minvar(gseE3.tlim(tint));
lmN = v(3,:);



ic = 4;
mspB = irf.tint('2015-10-16T10:33:20.000Z',0.2);
c_eval('Lmn = irf_norm(mean(dmpaB?brstRemOff.tlim(tint).data));',ic)

%%
tint = irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z');
[out,l,v] = irf_minvar(dmpaB1.tlim(tint));

L = v(1,:);
M = v(2,:); 
N = v(3,:);
%% N is fixed
N = lmN;
L = cross(cross(N,L),N);
M = cross(N,L);
[L;M;N]
%% L is fixed
L = Lmn;
N = cross(cross(L,lmN),L);
M = cross(N,L);
[L;M;N]
