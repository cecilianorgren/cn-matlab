% test yk.my_interf

%% load data to test on.
cd /Users/Cecilia/Research/Rongsheng/
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';
% Load internal burst probe data
IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
% Transform individual probe potentials to ISR2-system
diP1=cn.probe2isr2(IB(:,[1 2]),4,1);
diP2=cn.probe2isr2(IB(:,[1 3]),4,2);
diP3=cn.probe2isr2(IB(:,[1 4]),4,3);
diP4=cn.probe2isr2(IB(:,[1 5]),4,4);

%% try yk.my_interf on p3 and p4
% dont actually need to transform it to isr2 or something like that since
% it will only give a "?gonblicksbild"
[dt dtmax dtmin] = yk.my_interf(IB(:,[1 2]),IB(:,[1 3]),IB(:,[1 4]),IB(:,[1 5]));
