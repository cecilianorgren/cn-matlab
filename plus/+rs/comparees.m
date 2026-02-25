% Compares electric field obtained by different methods/ products.
% C4_CP_EFW_L1_IB, 
% C4_CP_EFW_L2_EB, 
% C4_CP_EFW_L3_E, 
% C4_CP_EFW_L1_P12/P34
%
%

% Go to directory and load time interval
cd /Users/Cecilia/Research/Rongsheng/
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';

%% 1. Load internal burst probe data
IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
% Transform individual probe potentials to ISR2-system
diP1=cn.probe2isr2(IB(:,[1 2]),4,1);
diP2=cn.probe2isr2(IB(:,[1 3]),4,2);
diP3=cn.probe2isr2(IB(:,[1 4]),4,3);
diP4=cn.probe2isr2(IB(:,[1 5]),4,4);
% Construct electric field
diEa = irf_add(-1,diP1,1,diP2); % (p2-p1)
diEb = irf_add(-1,diP3,1,diP4); % (p4-p3)
diE = irf_add(1,diEa,1,diEb); % total E
diEprobe = irf_tappl(diE,'*0.00212/0.044');

%% 1b. Load internal burst probe data
IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
% Transform individual probe potentials to ISR2-system
%diP1=cn.probe2isr2(IB(:,[1 2]),4,1);
%diP2=cn.probe2isr2(IB(:,[1 3]),4,2);
%diP3=cn.probe2isr2(IB(:,[1 4]),4,3);
%diP4=cn.probe2isr2(IB(:,[1 5]),4,4);
% Construct electric field
Ea = irf_add(-1,IB(:,[1 3]),1,IB(:,[1 2])); % (p2-p1)
Eb = irf_add(-1,IB(:,[1 5]),1,IB(:,[1 4])); % (p4-p3)
diEprobe12 = irf_tappl(Ea,'*0.00212/0.088');
diEprobe34 = irf_tappl(Eb,'*0.00212/0.088');

%% 2. Load internal burst electric field
diEcaa=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_EB','mat');


%% Compare the fields  
h=irf_plot({diEprobe,diEcaa},'comp');
for k=1:3; irf_legend(h(k),{'from probes','from caa'},[0.02 0.95]); end

%%
h=irf_plot({diEaprobe,diEbprobe,diEcaa},'comp');
for k=1:3; irf_legend(h(k),{'from probes 12','from probes 34','from caa'},[0.02 0.95]); end

%% Compare absolute value
Eabs1 = [diEcaa(:,1) sqrt(diEcaa(:,2).^2+diEcaa(:,3).^2)]; 
diEprobe = [diEaprobe diEbprobe(:,2)];
Eabs2 = [diEprobe(:,1) sqrt(diEprobe(:,2).^2+diEprobe(:,3).^2)];

g=irf_plot({Eabs1(:,[1 end]),Eabs2(:,[1 end])},'comp');
% Ok! the absolute value is the same, meaning the transformation to
% di-system is all thats left to fully construct E.

%% Transform to ISR2
diE12=cn.probe2isr2(diEaprobe,4,2);
diE34=cn.probe2isr2(diEbprobe,4,4);

%%
irf_plot({irf_add(1,diE12,1,diE34),diEcaa},'comp')
%% check if absoltue value is the same
Eabs3 = irf_add(1,irf_abs(diE12),1,irf_abs(diE34));
g=irf_plot({Eabs1(:,[1 end]),Eabs2(:,[1 end]),Eabs3(:,[1 end])},'comp');
%l=irf_plot(irf_add(1,Eabs2(:,[1 end]),-1,Eabs3(:,[1 end])));
% No, absolute value not conserved in transformation.
m=irf_plot(irf_multiply(1,Eabs2(:,[1 end]),1,Eabs3(:,[1 end]),-1));
