% mms data

% Define time int
tint = irf.tint('2015-05-15T00:40:00Z/2015-05-15T01:10:00Z');

%% Load sc position data
load /data/mms/irfu/mmsR.mat
epoTmp = EpochTT(R.time);
gsmR1 = [epoTmp.epochUnix R.gsmR1];

%% Load E, scP, B
% Takes a long time
for scId = 1
  fprintf('Loading MMS%d\n',scId);
  c_eval([...
    'E? = mms.db_get_ts(''mms?_edp_comm_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);'...
    'P? = mms.db_get_ts(''mms?_edp_comm_l2_scpot'',''mms?_edp_psp'',tint);'...
    'B? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_gsm_dmpa'',tint);'],...
    scId)
end
fprintf('Data loaded\n');


%% Try to make a wavelet

%%






