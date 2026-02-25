% line 600

%% [caa_is_get(100) : dsrc] ISDAT : empty return
probe = 32;
cl_id = 4;
[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);
            
%% [caa_is_get(100) : dsrc] ISDAT : empty return
probe = 12;
cl_id = 4;
[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);

%% [caa_is_get(98) : dsrc] ISDAT : no data
probe = 42;
cl_id = 3;
[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);
            
%%
probe = 32;
cl_id = 3;
[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);
            
%%
probe = 34;
cl_id = 3;
[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);