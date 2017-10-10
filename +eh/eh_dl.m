% Run probe-probe matching on eventual electron holes at the double layer.
cd /Users/Cecilia/Data/BM/20070831;
tint = toepoch([2007 08 31 10 17 00;2007 08 31 10 19 00])';
%%
c_eval('IB?=c_caa_var_get(''Data__C?_CP_EFW_L1_IB'',''mat'');',3:4);
c_eval('dobjIB?=c_caa_var_get(''Data__C?_CP_EFW_L1_IB'',''dobj'');',3:4);
% IB3 is the one at the double layer


%% try yk.my_interf on p3 and p4
% dont actually need to transform it toisr2 or something like that since
% it will only give a "?gonblicksbild"
IB=IB3;
[dt dtmax dtmin] = yk.my_interf(IB(:,[1 2]),IB(:,[1 3]),IB(:,[1 4]),IB(:,[1 5]));