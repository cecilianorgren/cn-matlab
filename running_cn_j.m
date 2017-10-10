% estimate j_par
ps3=c_caa_var_get('sc_pos_xyz_gse__C3_CP_FGM_FULL','mat');
ps4=c_caa_var_get('sc_pos_xyz_gse__C4_CP_FGM_FULL','mat');
c_eval('ps?=c_coord_trans(''gse'',''gsm'',ps?,''cl_id'',?);',3:4); % now in gsm
%bav=irf_add(-0.5,gsmB3,0.5,gsmB4);
%%
db=irf_add(-1,gsmB3,1,gsmB4);
dsc=irf_add(-1,ps3,1,ps4);
j=cn_j(db,dsc);
jpar = irf_dot(j,irf_norm(irf_add(0.5,gsmB3,0.5,gsmB4)));
figure(34);
jh=irf_plot({db,dsc,irf_tappl(j,'*1e9'),irf_tappl(jpar,'*1e9')}); 
tlim = toepoch([2007 08 31 10 17 30;2007 08 31 10 17 50])';
irf_zoom(jh,'x',tlim); irf_zoom(jh,'y')
irf_zoom(jh(3),'y',[-30 30])
tlim = toepoch([2007 08 31 10 17 38;2007 08 31 10 17 40])';
irf_pl_mark(jh,tlim)
ylabel(jh(1),'dB_{GSM} [nT]')
ylabel(jh(2),'dR_{GSM} [km]')
ylabel(jh(3),'j_{GSM} [nA/m^2]')
ylabel(jh(4),'j_{||,GSM} [nA/m^2]')
