%
%
% See also c_caa_distribution_data.m

% make spherical distribution functions
N=100;
[X,Y,Z] = sphere(N);
C=X(1:end-1,1:end-1);
C(:,:)=1;
surf(X,Y,Z,C)
colorbar

%% make own sphere
res=c_caa_construct_subspin_res_data('Data__C3_CP_PEA_3DXPH_DEFlux');
%%
c_eval('B?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3);

%%
r=10;
en_ind=25;
n_t=62;
ind_t=100:(100+n_t-1);
n_phi=n_t;
phi=torow(res.phiphi(1:n_phi));
theta=tocolumn([res.theta-15 res.theta(end)+15]);
ntheta=size(theta);

%theta=linspace(0,pi,ntheta);
%phi=linspace(0,2*pi,nphi);
X = r*sind(theta)*cosd(phi);
Y = r*sind(theta)*sind(phi);
Z = r*cosd(theta)*cosd(phi*0);
C = res.data(ind_t,:,en_ind);
surf(X,Y,Z,C'); hold(gca,'on');
colorbar

B=cn_toepoch(res.tt(1),B3);
B=irf_norm(B(2:end));
s=1.2;
plot3(s*r*B(1)*[-1 1],s*r*B(2)*[-1 1],s*r*B(3)*[-1 1])

axis square
axis equal
