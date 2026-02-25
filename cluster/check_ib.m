% Why is this not working?
cd /Users/Cecilia/Research/Rongsheng/
sc = 4;
probes = 1:4;
t1=[2003 08 17 16 50 00];
t2=[2003 08 17 17 00 00];
tint=toepoch([t1;t2])';

%%
c_eval('IB?=c_caa_var_get(''Data__C?_CP_EFW_L1_IB'',''mat'');',sc);

c_eval('P12_ibderv?=irf_add(-1,IB?(:,[1 2]),1,IB?(:,[1 3]));',sc);
c_eval('P34_ibderv?=irf_add(-1,IB?(:,[1 4]),1,IB?(:,[1 5]));',sc);

c_eval('E12_ibderv?=irf_tappl(P12_ibderv?,''/88'');',sc);
c_eval('E34_ibderv?=irf_tappl(P34_ibderv?,''/88'');',sc);

%% Without c_eval
if 0
    IB=c_caa_var_get('Data__C4_CP_EFW_L1_IB','mat');
    P12_ibderv=irf_add(-1,IB(:,[1 2]),1,IB(:,[1 3]));
    P34_ibderv=irf_add(-1,IB(:,[1 4]),1,IB(:,[1 5]));
    E12_ibderv=irf_tappl(P12_ibderv,''/88'');
    E34_ibderv=irf_tappl(P34_ibderv,''/88'');
end
%%
c_eval('[tt?,phase_data?] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], tint(1), diff(tint));',sc);
