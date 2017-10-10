% chinese guy
% Rongsheng Wang: see Questions on electron holes.ppt
%

cd /Users/Cecilia/Research/Rongsheng/
%% Download data to compare analysis
% Need to see what debye length is, therefore need temperature
%
t1=[2003 08 17 15 00 00];
t2=[2003 08 17 18 00 00];
tint=toepoch([t1;t2])';
%%
sc = 4;
probes = 1:4;
%c_eval('caa_download(tint,''C?_CP_PEA_MOMENTS'')',sc);
%%
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';

%c_eval('caa_download(tint,''C?_CP_EFW_L2_EB'');',sc);
%c_eval('caa_download(tint,''C?_CP_EFW_L1_P12'');',sc);
%c_eval('caa_download(tint,''C4_CP_EFW_L1_P?'');',1:4);
%c_eval('caa_download(tint,''C?_CP_EFW_L1_P34'');',sc);
%c_eval('caa_download(tint,''C?_CP_EFW_L1_IB'');',sc);
%c_eval('caa_download(tint,''C?_CP_EFW_L2_BB'');',sc);
%c_eval('caa_download(tint,''C?_CP_EFW_L3_SFIT'');',sc); 
%c_eval('caa_download(tint,''C?_CP_AUX_SPIN_TIME'');',sc); 
%c_eval('caa_download(tint,''C?_CH_AUX_SPIN_AXIS'');',sc); 
%% Check for internal burst data
%irf.tt_create_efw_internal_burst_list;

%% Load data
c_eval('[tt?,phase_data?] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], tint(1), diff(tint));',sc);
c_eval('phase?=[tt? phase_data?];',sc);
c_eval('spintime?=c_caa_var_get(''spin_period__C?_CP_AUX_SPIN_TIME'',''mat'');',sc);
c_eval('spinaxis?=c_caa_var_get(''__C?_CH_AUX_SPIN_AXIS'',''mat'');',sc);
c_eval('Te?MK=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sc);
c_eval('Te?=[Te?MK(:,1) irf_convert(Te?MK(:,2),''MK2eV'')];',sc)
c_eval('P12_?=c_caa_var_get(''P12__C?_CP_EFW_L1_P12'',''mat'');',sc);
c_eval('P34_?=c_caa_var_get(''P34__C?_CP_EFW_L1_P34'',''mat'');',sc);
c_eval('P?=c_caa_var_get(''P?__C4_CP_EFW_L1_P?'',''mat'');',1:4);
%c_eval('[caaVar,caaDobj,P34_?]=c_caa_var_get(''P34__C?_CP_EFW_L1_P34'',''mat'');',sc);
c_eval('Eib?=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sc);
c_eval('Bib?=c_caa_var_get(''B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB'',''mat'');',sc);
c_eval('diB?=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sc);
c_eval('IB?=c_caa_var_get(''Data__C?_CP_EFW_L1_IB'',''mat'');',sc);
c_eval('gseVsc?=c_caa_var_get(''sc_v_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',sc);
%% Change coordinate ssytem
c_eval('diVsc?=c_coord_trans(''gse'',''isr2'',gseVsc?,''cl_id'',?);',sc);
%% Make P12, P34
c_eval('P12_derv?=irf_add(-1,P1,1,P2);',sc);
c_eval('P34_derv?=irf_add(-1,P3,1,P4);',sc);
c_eval('P12_ibderv?=irf_add(-1,IB?(:,[1 2]),1,IB?(:,[1 3]));',sc);
c_eval('P34_ibderv?=irf_add(-1,IB?(:,[1 4]),1,IB?(:,[1 5]));',sc);
%% Check if potential difference derived from internal burst data is 
% compatible with calculated electric field
c_eval('E12_ibderv?=irf_tappl(P12_ibderv?,''/88'');',sc);
c_eval('E34_ibderv?=irf_tappl(P34_ibderv?,''/88'');',sc);


%% Convert data
% probe potential difference to electric field
% divide by 88 m = 0.088 km to get E in mV/m;
c_eval('E12_?=irf_tappl(P12_?,''/0.088'');',sc);
c_eval('E34_?=irf_tappl(P34_?,''/0.088'');',sc);

c_eval('E12_derv?=irf_tappl(P12_derv?,''/0.088'');',sc);
c_eval('E34_derv?=irf_tappl(P34_derv?,''/0.088'');',sc);
%% Derive phase between isr-coordinate system, B and probe pairs
probes=1:4;
% interpolated time series of phase
% goes from 0 to 360
c_eval('phase?=c_phase(diB?(:,1),[tt? phase_data?]);',sc); 
%%
c_eval('phase_p1=irf_tappl(phase?,''/180*pi + 3*pi/4'');',sc);
phase_p3 = phase_p1 - pi/2 ;
phase_p2 = phase_p1 + pi   ;
phase_p4 = phase_p1 + pi/2 ;
%        phase_heea=data.phase/180*pi-(30)/180*pi;
%        phase_leea=phase_heea+pi;
%        phase_rapid=data.phase/180*pi + 60.167/180*pi; % rapid phase
%        phase_sunsensor=data.phase/180*pi + 26.367/180*pi; % the location o fsun sensor
%% make plot of phases
h=irf_plot(6);
hca=irf_panel('p_sc4'); irf_plot(hca,phase4); ylabel(hca,'phase sc4');
hca=irf_panel('p1'); irf_plot(hca,phase_p1); ylabel(hca,'p1');
hca=irf_panel('p2'); irf_plot(hca,phase_p2); ylabel(hca,'p2');
hca=irf_panel('p3'); irf_plot(hca,phase_p3); ylabel(hca,'p3');
hca=irf_panel('p4'); irf_plot(hca,phase_p4); ylabel(hca,'p4');
hca=irf_panel('p1234'); irf_plot({phase_p1,phase_p2,phase_p3,phase_p4},'comp'); ylabel(hca,'all phases');
irf_zoom(h,'x',tint)
irf_plot_axis_align
%% time series of probe direction (going from origion (0,0,0)) in dsc 
% coordinate system, which is same as sr2
rp1=[phase_p1(:,1) 44*cos(phase_p1(:,2)) 44*sin(phase_p1(:,2)) phase_p1(:,2)*0]; % in DS reference frame
rp2=[phase_p2(:,1) 44*cos(phase_p2(:,2)) 44*sin(phase_p2(:,2)) phase_p1(:,2)*0];
rp3=[phase_p3(:,1) 44*cos(phase_p3(:,2)) 44*sin(phase_p3(:,2)) phase_p1(:,2)*0];
rp4=[phase_p4(:,1) 44*cos(phase_p4(:,2)) 44*sin(phase_p4(:,2)) phase_p1(:,2)*0];
%% transforming to isr2
c_eval('diRp?=c_coord_trans(''dsc'',''isr2'',rp?,''cl_id'',sc);',probes);
%% Try to make E_ISR2
c_eval('dip?=cn.probe2isr2(IB4(:,[1 (?+1) ]),4,?);',probes);
c_eval('dir?=cn.probe2isr2([IB4(:,1) ones(size(IB4))*44],4,?);',probes);

%% something wrong here, memory issues, dividing a matrix by another matrix,
% matrix-wise also creates a humongous matrix 90000x90000?
%E = (dip2-dip1)/(dir2-dir1);
dr=irf_add(-1,dir1,1,dir2);
dp=irf_add(-1,dip1,1,dip2);
E = irf_multiply(1,dp,1,dr,-1);

h=irf_plot({irf_abs(dr),irf_abs(dp),E});
yLabels={'r_{p2}-r_{p1}','p_{p2}-p_{p1}','(p_{p2}-p_{p1})/(r_{p2}-r_{p1})'};
for ii=1:3; ylabel(h(ii),yLabels{ii}); end
irf_plot_axis_align
%% Compare premade products p12 and p1-p2 / p2-p1
h=irf_plot({dp,P12_4});
yLabels={'p2-p1','P12'};
for ii=1:2; ylabel(h(ii),yLabels{ii}); end

%%
%t=get_time(1);
c_eval('E_corr?=irf_cross(diVsc?,irf_tappl(diB?,''*1e3*1e3*1e-9''));',sc);
%E_corr4=irf_cross(diVsc4,irf_tappl(diB4,'*1e-9');
irf_plot(E_corr4)
%% take angle with diB?
c_eval('Bdotp?=irf_dot(irf_norm(diRp?),irf_norm(diB4));',probes);
c_eval('angle_p?=[Bdotp?(:,1) acosd(Bdotp?)];',probes);
%% make plot of angles
h=irf_plot(2);
hca=irf_panel('1'); irf_plot(hca,irf_dot(irf_norm(diRp1),irf_norm(diB4))); ylabel(hca,'B dot p1');
hca=irf_panel('2'); irf_plot(hca,acosd(irf_dot(irf_norm(diRp1),irf_norm(diB4)))); ylabel(hca,'acosd( B dot p1)');
irf_zoom(h,'x',tint)
irf_plot_axis_align

%% make plot of angles
h=irf_plot(4);
hca=irf_panel('1'); irf_plot(hca,angle_p1); ylabel(hca,'B / p1');
hca=irf_panel('2'); irf_plot(hca,angle_p2); ylabel(hca,'B / p2');
hca=irf_panel('3'); irf_plot(hca,angle_p3); ylabel(hca,'B / p3');
hca=irf_panel('4'); irf_plot(hca,angle_p4); ylabel(hca,'B / p4');
irf_zoom(h,'x',tint)
irf_plot_axis_align

%% make angle between x_isr, y_isr2 and probepair 12, 34
x12=[diRp1(:,1) acosd(irf_dot(diRp1,[1 0 0]))]; % x 1
y12=[diRp1(:,1) acosd(irf_dot(diRp1,[0 1 0]))]; % x 1
x34=[diRp1(:,1) acosd(irf_dot(diRp3,[1 0 0]))]; % x 1
y34=[diRp1(:,1) acosd(irf_dot(diRp3,[0 1 0]))]; % x 1

xy12=[diRp1(:,1) acosd(irf_dot(irf_norm(diRp1),[1 0 0])) acosd(irf_dot(irf_norm(diRp1),[0 1 0]))]; % x 1
xy34=[diRp1(:,1) acosd(irf_dot(irf_norm(diRp3),[1 0 0])) acosd(irf_dot(irf_norm(diRp3),[0 1 0]))]; % x 1


%% Make some sort of plot
h=irf_plot(5);
hca=irf_panel('E'); irf_plot(hca,Eib4);  ylabel(hca,'E'); legend(hca,'x','y')
hca=irf_panel('P12'); irf_plot(hca,P12_4); ylabel(hca,'P12');
hca=irf_panel('E12'); irf_plot(hca,E12_4); ylabel(hca,'E12');
hca=irf_panel('P34'); irf_plot(hca,P34_4); ylabel(hca,'P34');
hca=irf_panel('E34'); irf_plot(hca,E34_4); ylabel(hca,'E34');
irf_zoom(h,'x',tint)
irf_plot_axis_align
%% Make compare plot of electric fields
h=irf_plot(6);
hca=irf_panel('E'); irf_plot(hca,Eib4); ylabel(hca,'E [mV/m]'); irf_legend(hca,{'E_{x,ISR2}','E_{y,ISR2}'},[0.02 0.9])
hca=irf_panel('E12/E'); irf_plot(hca,{E12_4,Eib4(:,[1 2])},'comp'); ylabel(hca,'E [mV/m]'); irf_legend(hca,{'P12/88m','E_{x,ISR2}'},[0.02 0.9])
hca=irf_panel('p12xy'); irf_plot(hca,xy12); ylabel(hca,'angle of booms 12 in isr2'); irf_legend(hca,{'angle to x_{ISR2}','angle to x_{ISR2}'},[0.02 0.9])
hca=irf_panel('E12/E'); irf_plot(hca,{E34_4,Eib4(:,[1 2])},'comp'); ylabel(hca,'E [mV/m]'); irf_legend(hca,{'P12/88m','E_{x,ISR2}'},[0.02 0.9])
hca=irf_panel('p34xy'); irf_plot(hca,xy34); ylabel(hca,'angle of booms 23 in isr2'); irf_legend(hca,{'angle to x_{ISR2}','angle to x_{ISR2}'},[0.02 0.9])
hca=irf_panel('phase'); irf_plot(hca,phase4); ylabel(hca,'Phase');
irf_zoom(h,'x',tint)
irf_plot_axis_align
%% Make compare plot of P12 / P1-P2 and P34 / P3-P4
h=irf_plot(2);
hca=irf_panel('P12'); irf_plot(hca,{P12_4,P12_derv4},'comp'); ylabel(hca,'P12 [mV/m]'); irf_legend(hca,{'P12','P2-P1'},[0.02 0.9])
hca=irf_panel('P34'); irf_plot(hca,{P34_4,P34_derv4},'comp'); ylabel(hca,'P12 [mV/m]'); irf_legend(hca,{'P12','P4-P3'},[0.02 0.9])



irf_zoom(h,'x',tint)
irf_plot_axis_align
