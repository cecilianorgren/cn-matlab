% emf_illustrations

%% problem 1.7, draw the potential
phi=@(ra)(exp(-ra)./ra);
phi_small=@(ra)(4*pi*dirac(ra));
ra=0:0.01:2;
plot(ra,phi(ra))
xlabel('r/a');ylabel('phi/phi_0')


%%
t1=[2007 08 31 10 18 00];t2=[2007 08 31 10 19 00];
b3=irf_tlim(gsmB3,[toepoch(t1) toepoch(t2)]);
b4=irf_tlim(gsmB4,[toepoch(t1) toepoch(t2)]);
e3=irf_tlim(gsmE3,[toepoch(t1) toepoch(t2)]);
e4=irf_tlim(gsmE4,[toepoch(t1) toepoch(t2)]);
d3=irf_tlim(gsmPos3,[toepoch(t1) toepoch(t2)]);
d4=irf_tlim(gsmPos4,[toepoch(t1) toepoch(t2)]);

B0=irf_norm(mean([b3(:,2:4); b4(:,2:4)]));
%
d0=mean(irf_add(1,d3,-1,d4));
d=d0(2:4); clear d0;
vec=[0 1 0];
%vec=irf_norm(d);
fb3=irf_newxyz(b3,irf_norm(cross(vec,B0)),0,B0);
fb4=irf_newxyz(b4,irf_norm(cross(vec,B0)),0,B0);
fe3=irf_newxyz(e3,irf_norm(cross(vec,B0)),0,B0);
fe4=irf_newxyz(e4,irf_norm(cross(vec,B0)),0,B0);
fd3=irf_newxyz(d3,irf_norm(cross(vec,B0)),0,B0);
fd4=irf_newxyz(d4,irf_norm(cross(vec,B0)),0,B0);

fj=cn_j(irf_add(1,fb3,-1,fb4),irf_add(1,fd3,-1,fd4));

%
h=irf_plot({fb3,fe3,fe4,fj})
irf_legend(h(1),{'x','y','z'},[0.98 0.95]);
irf_legend(h(2),{'x','y','z'},[0.98 0.95]);
irf_legend(h(3),{'x','y','z'},[0.98 0.95]);
irf_legend(h(4),{'x','y','z'},[0.98 0.95]);
ylabel(h(1),'fb3')
ylabel(h(2),'fe3')
ylabel(h(3),'fe4')
ylabel(h(4),'fj')