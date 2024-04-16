pd_orig = ePDist1;
b_orig = dmpaB1;
v_orig = dbcsVe1;
fact_orig = facTe1;
n_orig = ne1;


time = EpochTT('2017-07-03T05:26:42.679180175Z');
id = pd_orig.time.tlim(time + 0.5*0.03*[-1 1]);

id = id-250;
pd = pd_orig(id);
b = b_orig.resample(pd);
v = v_orig.resample(pd);
n = n_orig.resample(pd);
fact = fact_orig.resample(pd);
tpar = fact.xx;
tperp = (fact.yy + fact.zz)/2;



pd_omni = pd.omni;
pd_fpar = pd.reduce('1D',b,'vg',-6000:100:6000);

A_kappa = @(N,k,theta) N*pi*k*theta.^(-1/2)*gamma(k+1)/gamma(k-1/2);
fit_fun_kappa = @(v,N,k,theta) A_kappa(N,k,theta).*(1 + v.^2./k./theta.^2).^(-k-1);

A_flattop = @(N,k,theta) (1/2/pi)*N*theta^(-1)*gamma(1/k)/gamma(1+3/2/k)/abs(gamma(-1/2/k));
fit_fun_flattop = @(v,N,k,theta) A_flattop(N,k,theta).*(1 + (v./theta).^(2*k)).^((-k-1)/k);
  


%fit = ph_fit_kappa(pd_fpar,'x0',[0.4*1e6 2 2e5],'plot',1,'weight',[1 1 1]);
fit = ph_fit_kappa(pd_fpar,'x0',[10e9 3 4e6],'plot',1,'weight',[1 1 1]);


nrows = 4;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,pd_fpar.depend{1}*1e-3,pd_fpar.data,pd_fpar.depend{1}*1e-3,fit.f)


hca = h(isub); isub = isub + 1;
plot(hca,pd_fpar.depend{1}*1e-3,fit.f)


hca = h(isub); isub = isub + 1;
plotyy(hca,pd_fpar.depend{1}*1e-3,fit_fun_kappa(pd_fpar.depend{1}*1e3,1,1,2e6),...
         pd_fpar.depend{1}*1e-3,fit_fun_flattop(pd_fpar.depend{1}*1e3,1,2,2e6))

hca = h(isub); isub = isub + 1;
plot(hca,pd_fpar.depend{1}*1e-3,fit_fun_flattop(pd_fpar.depend{1},fit.X(1),fit.X(2),fit.X(3)))