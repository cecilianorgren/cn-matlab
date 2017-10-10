% Go to directory and load time interval
cd /Users/Cecilia/Research/Rongsheng/
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';

%% Load internal burst probe data
IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
% Construct electric field
prelE12 = irf_add(-1,IB(:,[1 3]),1,IB(:,[1 2])); % "(p2-p1)"
prelE34 = irf_add(-1,IB(:,[1 5]),1,IB(:,[1 4])); % "(p4-p3)"
E12 = irf_tappl(prelE12,'*0.00212/0.088'); % transform to right unities
E34 = irf_tappl(prelE34,'*0.00212/0.088'); % transform to right unities

%% Transform to ISR2 
% (Do this in order to compare with CAA product. If you just look at 
% individual holes you just need to know the instantaneous value)
diE12 = cn.probe2isr2(E12,4,2);
diE34 = cn.probe2isr2(E34,4,4);
% Add the contributions from each probe pair
diE1234 = irf_add(1,diE12,1,diE34);

%% Load CAA internal burst electric field
diEcaa=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_EB','mat');

%% Compare the fields  
h=irf_plot({diE1234,diEcaa},'comp');
for k=1:3; irf_legend(h(k),{'from probes','from caa'},[0.02 0.95]); end