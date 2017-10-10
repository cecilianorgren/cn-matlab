%% Go to correct directory and load time interval
cd /Users/Cecilia/Research/Rongsheng/
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';
%% Load internal burst probe data
IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
%% Transform individual probe potentials to ISR2-system
diP1=cn.probe2isr2(IB(:,[1 2]),4,1);
diP2=cn.probe2isr2(IB(:,[1 3]),4,2);
diP3=cn.probe2isr2(IB(:,[1 4]),4,3);
diP4=cn.probe2isr2(IB(:,[1 5]),4,4);
%% Transform individual probe coordinates to ISR2-system
diR1=cn.probe2isr2([IB(:,1), ones(size(IB(:,1)))*44],4,1);
diR2=cn.probe2isr2([IB(:,1), ones(size(IB(:,1)))*44],4,2);
diR3=cn.probe2isr2([IB(:,1), ones(size(IB(:,1)))*44],4,3);
diR4=cn.probe2isr2([IB(:,1), ones(size(IB(:,1)))*44],4,4);
%% Take difference in probe coordinates.
didP12=irf_add(-1,diP1,1,diP2);
didP34=irf_add(-1,diP3,1,diP4);
%% Take difference in probe potentials.
didR12=irf_add(-1,diR1,1,diR2);
didR34=irf_add(-1,diR3,1,diR4);
%% Calculate electric field dy dividing potential difference by
% spatial separation distance.
diE12 = irf_multiply(1,didP12,1,didR12,-1);
diE34 = irf_multiply(1,didP34,1,didR34,-1);
%% Add the electric field from thw two probe pairs to obtain the
% total electric field.
diE = irf_add(1,diE12,1,diE34);
%% Compare with already calculated version
%diEib=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_EB','mat');
irf_plot({diE,diEib},'comp')
% Results are that they are not at all the same

