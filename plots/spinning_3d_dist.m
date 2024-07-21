mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');
%db_info = datastore('mms_db');

%% Torbert event 
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

%[enflux_new, enflux_BG, idist_new, idist_BG, Ni_new, gseVi_new, gsePi_new, ...
%  Ni_bg, EnergySpectr_bg, Pres_bg, EnergySpectr_bg_self]= mms.remove_ion_penetrating_radiation_bg(iPDist3);

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

%%

L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn_edr = [L;M;N];

L_vi = -[-0.8906    0.4548    0.0045];
M_vi = [ 0.4539    0.8893   -0.0559];
N_vi = -[-0.0294   -0.0477   -0.9984];
lmn_vi = [L_vi; M_vi; N_vi];


L_gse = [1 0 0];
M_gse = [0 1 0];
N_gse = [0 0 1];
lmn_gse = [L_gse; M_gse; N_gse];

L_gse = [1 0 0];
M_gse = [0 1 -0.2];
%N_gse = [0 0 1];
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

lmn = lmn_gse;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);


c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''E LMN'';',ic)
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn_vi''; mvaVi?.name = ''Vi LMN'';',ic)

c_eval('tsLgse? = irf.ts_vec_xyz(iPDist?.time,repmat(L,iPDist?.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(iPDist?.time,repmat(M,iPDist?.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(iPDist?.time,repmat(N,iPDist?.length,1));',ic)

c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt?,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt?,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt?,-1);',ic)

%% Plot
c_eval('pdist_movmean = iPDist?.movmean(11,''RemoveOneCounts'',iPDist?_counts);',ic)

dt = 5;
elim = [200 Inf];
time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
time = time + dt;
tint_dist = time + 1*0.5*0.150*[-1 1];

pdist = pdist_movmean.tlim(tint_dist);


t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
%c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
%Ldsl = Tdsl(1,:);
%Mdsl = Tdsl(2,:);
%Ndsl = Tdsl(3,:);

%lmn = [Ldsl, Mdsl, Ndsl];

lmn = [1 0 0; 0 1 0; 0 0 1];

h = subplot(1,1,1);
isub = 1;

nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
iso_values = 20e-28;
vlim = 3000;
if 1 % isuorface
  hca = h(isub); isub = isub + 1;
  hca.ColorOrder = pic_colors('matlab');
  hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
  %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
  %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
  c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
  axis(hca,'square')
  
  hca.XLim = vlim*[-1 1];
  hca.YLim = vlim*[-1 1];
  hca.ZLim = vlim*[-1 1];
 
  camlight(gca,0,0)
  
  %h(isub-4).Title = hca.Title;
  %h(isub-4).Title.FontSize = 8;
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_M (km/s)';
  hca.ZLabel.String = 'v_N (km/s)';
  hca.Title = [];
end

%% Annotate
% arrow for tip of arrow head/triangle
xx = [500 500];
yy = [1500 1000];
zz = [1200 1200];
[harr,xxx,yyy] = arrow([xx(1) yy(1) zz(1)],[xx(2) yy(2) zz(2)],'color','k');

%% Make movie
vidObj = VideoWriter([printpath ,'3D_turning_4h'],'MPEG-4');
open(vidObj)


nAngles = 36;
%angles = linspace(0,360,nAngles);

set(gcf,'color','white');
angles = 0:2:358;
nAngles = numel(angles);
for ia = 1:nAngles
  viewAxis1 = [0 0 1]; 
  
  viewAxis1 = viewAxis1/norm(viewAxis1);
  viewAxis3 = cross(viewAxis1,[1 0 0]); viewAxis3 = viewAxis3/norm(viewAxis3);
  viewAxis2 = cross(viewAxis3,viewAxis1); viewAxis2 = viewAxis2/norm(viewAxis2);
  viewDir = viewAxis2*cosd(angles(ia)) + viewAxis3*sind(angles(ia));
  %camDir = viewAxis2*cosd(angles(ia)) + viewAxis3*sind(angles(ia));
  view(viewDir)
  delete(findobj(gcf,'type','light'))
  
  camlight(gca,0,0)

  hca.Visible = 'off';
  %hca.OuterPosition = [0.1 0.1 0.8 0.8];
  %hca.PositionConstraint = 'innerposition';
  %hca.PositionConstraint = 'outerposition';
  %hca.Projection = "perspective";
  %hca.InnerPosition = [0.2170    0.1990    0.6975    0.7335] - [0.05 0.05 -0.05 0.05];
  %hca.InnerPosition = [0 0 1 1];
  %hca.InnerPosition = [0.1 0.1 0.8 0.8];
  %hca.InnerPosition = [0.1 0.1 0.8 0.8];
  
  axis(hca,'square')
  axis(hca,'equal')
  camp = campos;
  %norm(camp)
  hca.CameraViewAngle = 12;
  hca.CameraViewAngleMode = 'manual';
  pause(0.1)

  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end

close(vidObj)











