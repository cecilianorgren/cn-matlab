ic = 1:4;
if 0
db_info = datastore('mms_db');  
dbroot = db_info.local_file_db_root;


interval = '20151016103254';
yyyy = interval(1:4);
mm = interval(5:6);
dd = interval(7:8);
starttime = interval(9:end);

tmp = ['tmpDataObj? = dataobj(''' dbroot '/mms?/fpi/brst/l2/des-moms/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-moms_' interval '_v3.1.0.cdf'');'];
c_eval(tmp,ic);
c_eval('ne? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_numberdensity_brst''));',ic);

tint = irf.tint(ne1.time.start.utc,ne1.time.stop.utc);
tintl = tint+[-60 60];
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

R  = mms.get_data('R_gse',tintl);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',ic);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',ic);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);
Bmagav = Bxyzav.abs;

bav = Bxyzav/Bxyzav.abs;
end
%%
c_eval('Rxyz? = mvaR?.resample(mvaB1);',1:4)
c_eval('Bxyz? = mvaB?.resample(mvaB1);',1:4)
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);
Bmagav = Bxyzav.abs;

bav = Bxyzav/Bxyzav.abs;

[curvb1,~]=c_4_grad('Rxyz?','Bxyz?','curvature');

c_eval('b? = Bxyz?/Bxyz?.abs;',ic);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);

[gradb,~]=c_4_grad('Rxyz?','b?','grad');
curvbx2 = bav.x.data.*gradb.data(:,1)+bav.y.data.*gradb.data(:,4)+bav.z.data.*gradb.data(:,7);
curvby2 = bav.x.data.*gradb.data(:,2)+bav.y.data.*gradb.data(:,5)+bav.z.data.*gradb.data(:,8);
curvbz2 = bav.x.data.*gradb.data(:,3)+bav.y.data.*gradb.data(:,6)+bav.z.data.*gradb.data(:,9);
curvb2 = irf.ts_vec_xyz(Bxyz1.time,[curvbx2 curvby2 curvbz2]);

%% 

[gradB,~]=c_4_grad('Rxyz?','Bxyz?','grad');

A = Bxyzav.x.data.*Bxyzav.x.data.*gradB.data(:,1);
B = Bxyzav.x.data.*Bxyzav.y.data.*gradB.data(:,2);
C = Bxyzav.x.data.*Bxyzav.z.data.*gradB.data(:,3);
D = Bxyzav.y.data.*Bxyzav.x.data.*gradB.data(:,4);
E = Bxyzav.y.data.*Bxyzav.y.data.*gradB.data(:,5);
F = Bxyzav.y.data.*Bxyzav.z.data.*gradB.data(:,6);
G = Bxyzav.z.data.*Bxyzav.x.data.*gradB.data(:,7);
H = Bxyzav.z.data.*Bxyzav.y.data.*gradB.data(:,8);
I = Bxyzav.z.data.*Bxyzav.z.data.*gradB.data(:,9);
correction = A+B+C+D+E+F+G+H+I;

curvbx3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,1)+Bxyzav.y.data.*gradB.data(:,4)+Bxyzav.z.data.*gradB.data(:,7));
curvby3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,2)+Bxyzav.y.data.*gradB.data(:,5)+Bxyzav.z.data.*gradB.data(:,8));
curvbz3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,3)+Bxyzav.y.data.*gradB.data(:,6)+Bxyzav.z.data.*gradB.data(:,9));
curvbx3 = curvbx3 - (1./Bmagav.data).^4.*(Bxyzav.x.data.*correction);
curvby3 = curvby3 - (1./Bmagav.data).^4.*(Bxyzav.y.data.*correction);
curvbz3 = curvbz3 - (1./Bmagav.data).^4.*(Bxyzav.z.data.*correction);
curvb3 = irf.ts_vec_xyz(Bxyz1.time,[curvbx3 curvby3 curvbz3]);

RC = 1/curvb3.abs;
RCE = RC/6400;


h = irf_plot({Bxyzav,curvb1,curvb2,curvb3,RC});
set(h(5),'yscale','log');
%irf_plot({Bxyzav,dot(curvb1,bav),dot(curvb2,bav),dot(curvb3,bav)})