ic = 1:4;

db_info = datastore('mms_db');  
dbroot = db_info.local_file_db_root;


interval = '20151112071854';
yyyy = interval(1:4);
mm = interval(5:6);
dd = interval(7:8);
starttime = interval(9:end);

%tmp = ['tmpDataObj? = dataobj(''' dbroot '/mms?/fpi/brst/l2/des-moms/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-moms_' interval '_v3.1.0.cdf'');'];
tmp = ['tmpDataObj? = dataobj(''' dbroot '/mms?/fpi/brst/l2/des-moms/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-moms_' interval '_v3.1.0.cdf'');'];
c_eval(tmp,ic);
c_eval('ne? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_numberdensity_brst''));',ic);
%ls
tint = irf.tint(ne1.time.start.utc,ne1.time.stop.utc);
tintl = tint+[-60 60];
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

R  = mms.get_data('R_gse',tintl);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',ic);
c_eval('bxyz? = Bxyz?/Bxyz?.abs;',ic);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',ic);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);
Bmagav = Bxyzav.abs;

bav = Bxyzav/Bxyzav.abs;

%% IRFU-MATLAB
[curvb1,~]=c_4_grad('Rxyz?','Bxyz?','curvature'); % (b*grad)b
curvb1.name = 'c 4 grad';
[curvb_1,~]=c_4_grad('Rxyz?','bxyz?','curvature'); % (b*grad)b
curvb_1.name = 'c 4 grad';


% Make units vectors before calling c_4_grad, multiply the grad-terms with
% the average vector.
c_eval('b? = Bxyz?/Bxyz?.abs;',ic);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);
[gradb,~]=c_4_grad('Rxyz?','b?','grad');
curvbx2 = bav.x.data.*gradb.data(:,1)+bav.y.data.*gradb.data(:,4)+bav.z.data.*gradb.data(:,7);
curvby2 = bav.x.data.*gradb.data(:,2)+bav.y.data.*gradb.data(:,5)+bav.z.data.*gradb.data(:,8);
curvbz2 = bav.x.data.*gradb.data(:,3)+bav.y.data.*gradb.data(:,6)+bav.z.data.*gradb.data(:,9);
curvb2 = irf.ts_vec_xyz(Bxyz1.time,[curvbx2 curvby2 curvbz2]);
curvb2.name = 'irfu-m (corr)';

%% Shen 2003, with correction

[gradB,~]=c_4_grad('Rxyz?','Bxyz?','grad'); % grad B

if 1
  dxBx = gradB.data(:,1);
  dxBy = gradB.data(:,2);
  dxBz = gradB.data(:,3);
  dyBx = gradB.data(:,4);
  dyBy = gradB.data(:,5);
  dyBz = gradB.data(:,6);
  dzBx = gradB.data(:,7);
  dzBy = gradB.data(:,8);
  dzBz = gradB.data(:,9);
  Bx = Bxyzav.x.data;
  By = Bxyzav.y.data;
  Bz = Bxyzav.z.data;
  B = Bxyzav.abs.data;
  
  Cxx = Bx.*Bx.*dxBx;
  Cxy = Bx.*By.*dxBy;
  Cxz = Bx.*Bz.*dxBz;
  Cyx = By.*Bx.*dyBx;
  Cyy = By.*By.*dyBy;
  Cyz = By.*Bz.*dyBz;
  Czx = Bz.*Bx.*dzBx;
  Czy = Bz.*By.*dzBy;
  Czz = Bz.*Bz.*dzBz;
  correction = Cxx+Cxy+Cxz+Cyx+Cyy+Cyz+Czx+Czy+Czz;
  
  curv_correction_x = - (1./B).^4.*(Bx.*correction);
  curv_correction_y = - (1./B).^4.*(By.*correction);
  curv_correction_z = - (1./B).^4.*(Bz.*correction);
  
  curv_correction = irf.ts_vec_xyz(Bxyz1.time,[curv_correction_x curv_correction_y curv_correction_z]);
  curv_correction.name = 'Correction';
  
  if 0
    curvbx3 = (1./B).^2.*(Bx.*dxBx + By.*dyBx + Bz.*dzBx);
    curvby3 = (1./B).^2.*(Bx.*dxBy + By.*dyBy + Bz.*dzBy);
    curvbz3 = (1./B).^2.*(Bx.*dxBz + By.*dyBz + Bz.*dzBz);
  else   
    c_eval('curvbx3_? = (1./Bxyz?.abs.data).^2.*(Bxyz?.x.data.*dxBx + Bxyz?.y.data.*dyBx + Bxyz?.z.data.*dzBx);',1:4)
    c_eval('curvby3_? = (1./Bxyz?.abs.data).^2.*(Bxyz?.x.data.*dxBy + Bxyz?.y.data.*dyBy + Bxyz?.z.data.*dzBy);',1:4)
    c_eval('curvbz3_? = (1./Bxyz?.abs.data).^2.*(Bxyz?.x.data.*dxBz + Bxyz?.y.data.*dyBz + Bxyz?.z.data.*dzBz);',1:4)
    curvbx3 = (curvbx3_1 + curvbx3_2 + curvbx3_3 + curvbx3_4)/4;
    curvby3 = (curvby3_1 + curvby3_2 + curvby3_3 + curvby3_4)/4;
    curvbz3 = (curvbz3_1 + curvbz3_2 + curvbz3_3 + curvbz3_4)/4;
  end
  curvb3_nocorrection = irf.ts_vec_xyz(Bxyz1.time,[curvbx3 curvby3 curvbz3]); curvb3_nocorrection.name = 'no corr';

  curvbx3 = curvbx3 + curv_correction_x;
  curvby3 = curvby3 + curv_correction_y;
  curvbz3 = curvbz3 + curv_correction_z;
else  
  % correction to keep curvature perpendicular to average magnetic field:
  % B^(-4)*Bj*Bi*Bl*(di/Bl)
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
  curv_correction_x = - (1./Bmagav.data).^4.*(Bxyzav.x.data.*correction);
  curv_correction_y = - (1./Bmagav.data).^4.*(Bxyzav.y.data.*correction);
  curv_correction_z = - (1./Bmagav.data).^4.*(Bxyzav.z.data.*correction);
  curv_correction = irf.ts_vec_xyz(Bxyz1.time,[curv_correction_x curv_correction_y curv_correction_z]);
  curv_correction.name = 'Correction';
  
  curvbx3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,1)+Bxyzav.y.data.*gradB.data(:,4)+Bxyzav.z.data.*gradB.data(:,7));
  curvby3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,2)+Bxyzav.y.data.*gradB.data(:,5)+Bxyzav.z.data.*gradB.data(:,8));
  curvbz3 = (1./Bmagav.data).^2.*(Bxyzav.x.data.*gradB.data(:,3)+Bxyzav.y.data.*gradB.data(:,6)+Bxyzav.z.data.*gradB.data(:,9));
  curvb3_nocorrection = irf.ts_vec_xyz(Bxyz1.time,[curvbx3 curvby3 curvbz3]); curvb3_nocorrection.name = 'no corr';

  curvbx3 = curvbx3 - (1./Bmagav.data).^4.*(Bxyzav.x.data.*correction);
  curvby3 = curvby3 - (1./Bmagav.data).^4.*(Bxyzav.y.data.*correction);
  curvbz3 = curvbz3 - (1./Bmagav.data).^4.*(Bxyzav.z.data.*correction);
end



curvb3 = irf.ts_vec_xyz(Bxyz1.time,[curvbx3 curvby3 curvbz3]); 
curvb3.name = 'Curv B (Shen)';

RC = 1/curvb3.abs; RC.name = 'R_C (Shen)';
RCE = RC/6400;

%% Plot
tint = irf.tint('2015-11-12T07:19:20.50Z/2015-11-12T07:19:22.00Z');
h = irf_plot({Bxyzav,curvb1,curvb_1,curvb2,curvb3,curvb3_nocorrection,curv_correction,RC});
set(h(end),'yscale','log','ytick',10.^[-1:3],'ylim',[1e1 3e4]);
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
%irf_plot({Bxyzav,dot(curvb1,bav),dot(curvb2,bav),dot(curvb3,bav)})