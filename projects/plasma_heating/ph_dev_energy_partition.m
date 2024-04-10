function varargout = ph_dev_energy_partition(pd_orig,b_orig,v_orig,fact_orig)
% PH_DEV_ENERGY_PARTITION Quantifies energy partition of PDist
%
%   results = PH_DEV_ENERGY_PARTITION(PDist,dmpaB,dslV,facT);
%
%
%

%%

pd = pd_orig;
b = b_orig.resample(pd);
v = v_orig.resample(pd);
fact = fact_orig.resample(pd);
tpar = fact.xx;
tperp = (fact.yy + fact.zz)/2;


pd_omni = pd.omni;
pd_fpar = pd.reduce;


f_A = @(N,k,theta) N*pi*k*theta^(-3/2)*gamma(k+1)/gamma(k-1/2);
f_kappa = @(A,v,k,theta) f_A(N,k,theta)*(1 + v.^2/k/theta^2)^(-k-1);

varargout{1} = [];