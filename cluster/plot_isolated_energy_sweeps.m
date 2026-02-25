% plot only isolated energy sweeps
t1=[2007 08 31 10 18 42.00];
t2=[2007 08 31 10 18 43.00];
tint=[toepoch(t1) toepoch(t2)];
energy_sweeps=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_3DXH_PSD','mat');
res=c_caa_construct_subspin_res_data('Data__C3_CP_PEA_PITCH_3DXH_PSD');
energy_sweeps_dl=irf_tlim(energy_sweeps,tint);

sweep_mode=c_caa_var_get('Mode_SweepMode__C3_CP_PEA_PITCH_3DXH_PSD','mat');

[variable,dataobject,peace,dataunits]=c_caa_var_get('Data__C3_CP_PEA_PITCH_3DXH_PSD');

%% ok seems ok, or
% during one time step, which is one azimuthal angle bin, one energy sweep
% is made, which could encompass a number of polar angles.
for k=1:12
    loglog(squeeze(res.en),squeeze(res.data(10000,k,:)),'color',0.8*[rand,rand,rand]); hold on;
end
%%
for k=1:10
h(k)=subplot(2,5,k);
c_caa_distribution_function(h(k),'tint',[tint(1)+k*0.125 tint(1)+(k-1)*0.14],'C4_CP_PEA_PITCH_3DXH_PSD','polar')
end

%% New plot with more pitch angles
tstart = [2007 08 31 10 17 30.00];
tstop  = [2007 08 31 10 18 15.00];
dt     = 0;

for k=1:10
h(k)=subplot(2,5,k);
c_caa_distribution_function(h(k),'tint',tint(1),'C4_CP_PEA_PITCH_3DXH_PSD','polar')
end
%tint = [toepoch(t1) toepoch(t2)];

