% Make N bi-maxwellian fit to FPI distribution
%% Load and prepare data
% Spacecraft id
ic = 1;
units = irf_units;
T = 0.03; % electron, 0.15 for ions
% Time interval of event
tint_burst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_burst = tint_burst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file

% Time interval for figure
tint_dist = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint_dist = irf_time('2017-07-06T13:54:05.56Z','utc>EpochTT') + 0.5*T*[-1 1];

% Set up database
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS'); 
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

% Load data
% Necessary data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint_burst);',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint_burst+[20 0]));',ic)
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint_burst);',ic);

% Optional data, to specify starting guess
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint_burst,?);',ic)
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint_burst,?);',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint_burst,?);',ic);

%% Rebin distribution function into cartesian grid
eint = [000 40000];
vint = [-Inf Inf];
vg = (-50:1:50)*1e3;

c_eval('eDist = ePDist?.tlim(tint_dist);',ic)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist(1)).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;
%lmn = [ePara.data; ePerp1.data; ePerp2.data];

%orient = [ePara.data; ePerp1.data; ePerp2.data];
orient = [ePerp1.data; ePerp2.data; ePara.data];
lowerelim = 40;
nMC = 100e0;

% For 3D/2D/1D fit (can further be reduced to 1D or 2D) (quite slow, buggy?)
%ef3D_perp1perp2par = eDist.elim([40 inf]).rebin('cart',{vg,vg,vg},orient); % Rebins skymap into 3D cartesian grid

% For 2D/1D fit (can further be reduced to 1D)
ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);
ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);
ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);

% For 1D fit
ef1D_par = eDist.reduce('1D',ePara.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B
ef1D_perp1 = eDist.reduce('1D',ePerp1.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B
ef1D_perp2 = eDist.reduce('1D',ePerp2.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B

% Get some values for initial guess of search function
c_eval('n0 = ne?.resample(eDist);' ,ic)
c_eval('v0 = gseVe?*orient''; v0 = v0.resample(eDist);' ,ic)
c_eval('T0 = mms.rotate_tensor(gseTe?,''rot'',orient(1,:),orient(2,:),orient(3,:)); T0 = T0.resample(eDist);',ic) % xx component is par

%% Make fit, 1D
nMax = 1; % how many Maxwellians to fit to
%dv = (ePDistCart.depend{1}(2) - ePDistCart.depend{1}(1));
%f_obs = squeeze(sum(sum(ePDistCart.data,4),3))*dv;
%v = ePDistCart.depend{1}*1e3; % m/s
ef = ef1D_par;
f_obs = ef.data;
v = ef.depend{1}*1e3; % m/s
params0 = [n0.data v0.data T0.xx.data T0.yy.data T0.zz.data];
T0scalar = T0.trace.data/3; % eV
vt0 = sqrt(T0scalar*units.eV*2/units.me)/1000; % km/s

% 1D
params_obs = [n0.data*1e6 v0.z.data*1e3 vt0*1e3]; % m^-3, m/s, m/s
nPop = 1; params0 = params_obs;
nPop = 2; params0 = [params_obs params_obs]; %params0 = params0.*rand(1,numel(params0));
%nPop = 2; params0 = [0.03e6 -5000e3 5000e3 0.01e6 1000e3 15000e3];
%nPop = 3; params0 = repmat(params_obs,1,nPop); params0 = params0.*(0.5+0.5*rand(1,numel(params0)));

% fminsearch uses the Nelder-Mead simplex (direct search) method.
% I'm not sure if the large difference between parameters affects the 
% method. For example, vt/n = 2e2. If the method varies each parameter 
% which the same "step", vt will not be varied much...
options = optimset('TolFun',0.2);
doPlot = 0;
cost_function = @(params) costfunction_maxwellian(params,{v},f_obs,nPop,0);
[params,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,params0,options);

% Doing 3 populations form the start works poorly for some reason, so one
% can first try two, then use the results as input to second search
nPop = 3; params0 = [params params_obs];
cost_function = @(params) costfunction_maxwellian(params,{v},f_obs,nPop,1);
params = fminsearch(cost_function,params0,options);

%% Make fit, 2D
nMax = 1; % how many Maxwellians to fit to
%dv = (ePDistCart.depend{1}(2) - ePDistCart.depend{1}(1));
%f_obs = squeeze(sum(sum(ePDistCart.data,4),3))*dv;
%v = ePDistCart.depend{1}*1e3; % m/s
ef = ef2D_parperp1;
f_obs = squeeze(ef.data);
v = ef.depend{1}*1e3; % m/s

T0scalar = T0.trace.data/3; % eV
T0par = T0.zz.data; % eV
T0perp = (T0.xx.data+T0.yy.data)/2; % eV
vt0 = sqrt(T0scalar*units.eV*2/units.me)/1000; % km/s
vtpar = sqrt(T0par*units.eV*2/units.me)/1000; % km/s
vtperp = sqrt(T0perp*units.eV*2/units.me)/1000; % km/s
vtx0 = sqrt(T0.xx.data*units.eV*2/units.me)/1000; % km/s
vty0 = sqrt(T0.yy.data*units.eV*2/units.me)/1000; % km/s
vtz0 = sqrt(T0.zz.data*units.eV*2/units.me)/1000; % km/s
vdpar0  = v0.z.data;
vdperp0  = v0.x.data;


% 2D
params_obs = [n0.data*1e6/2 vdpar0 vdperp0 [vtpar vtperp]*1e3]; % m^-3, m/s, m/s
nPop = 1; params0 = params_obs;
nPop = 2; params0 = [params_obs params_obs]; %params0 = params0.*rand(1,numel(params0));


% fminsearch uses the Nelder-Mead simplex (direct search) method.
% I'm not sure if the large difference between parameters affects the 
% method. For example, vt/n = 2e2. If the method varies each parameter 
% which the same "step", vt will not be varied much...

% Save value of costfunction.
history = [];
history_param = [];

options = optimset('OutputFcn', @myoutput);
options = optimset([]);

% Option to plot output. 2-D plots for each time step is consuming quite a
% lot of resources and is only good or initial checks.
doPlot = 0;
cost_function = @(params) costfunction_maxwellian(params,{v,v},f_obs,nPop,doPlot);
[params,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,params0);

cost_function = @(params) costfunction_maxwellian(params,{v,v},f_obs,nPop,0);
[params,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,params);

% % Doing 3 populations form the start works poorly for some reason, so one
% % can first try two, then use the results as input to second search
% nPop = 3; params0 = [params params_obs];
cost_function = @(params) costfunction_maxwellian(params,{v,v},f_obs,nPop,1);
params = fminsearch(cost_function,params,options);


%% Try with constrained search 
% A is matrix, b is a vector
% A*params <= b
% keep n >0
% A*n <= b
% -1*n <= 0   <- should work
% keep vt,vd < c
% 1*vd <= c
%
% [1 0 0;
%  0 1 0;
%  0 0 1]
c = units.c;
A = [-1 1 1]; A = repmat(A,1,nPop);
b = [0 c c]; b = repmat(b,1,nPop);
A = eye(9);
A(1,1) = -1;
A(4,4) = -1;
A(7,7) = -1;
b = [0 c c 0 c c 0 c c];
lb = zeros(size(params));
ub = [200 c c 200 c c 200];

options = optimset('TolFun',1e-3);

%params = fmincon(cost_function,params0,A,b');
params = fmincon(cost_function,params0,[],[],[],[],lb,ub);


%cost_function = @(params) eh_costfunction(params,xdata,ydata,zdata,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez);

%% Auxiliary/help functions
function stop = plotoutput(x,optimvalues,state,plotiter)
% Option to plot every nth iteration, to save time, but still get an idea
% of progress
  stop = false;
  if isequal(state,'iter')
    if state == plotiter
      doPlot = 1;  
    else
      doPlot = 0;
    end
  end
end
function stop = myoutput(x,optimvalues,state)
  stop = false;
  if isequal(state,'iter')
    history = [history, x];
  end
end