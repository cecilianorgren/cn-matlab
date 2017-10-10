cd /Users/Cecilia/Data/20030929

t1=[2003 09 29 10 18 00];
t2=[2003 09 29 10 21 00];
tint=[toepoch(t1) toepoch(t2)];

%% Load data: E, B
% E: C1-C3 for ISR2, C4 for GSE
% B: C1-C4
c_eval('diB?=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');');
c_eval('diPos?=c_caa_var_get(''sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');');
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',1:3);
c_eval('gseE?=c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''mat'');',4);
c_eval('diE?=c_coord_trans(''gse'',''dsi'',gseE?,''cl_id'',?);',4);

%% Construct complete E, assuming ExB=0
c_eval('[diE?t,angle?t]=irf_edb(diE?,diB?,1);');
c_eval('[diE?p,angle?p]=irf_edb(diE?,diB?,89,''Epar'');');
%% Plot fields in ISR2 system
h=irf_plot(8);
isub=1;

if 1 % diB_X
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diB?',2);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diB_X')
end
if 1 % diB_Y
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diB?',3);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diB_Y')
end
if 1 % diB_Z
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diB?',4);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diB_Z')
end
if 1 % diE_X
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diE?',2);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diE_X')
end
if 1 % diE_Y
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diE?',3);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diE_Y')
end
if 1 % diE_Y
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diE?p',4);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diE_Z ExB=0')
end
if 1 % diE_Y
    hca=h(isub); isub=isub+1;
    c_pl_tx(hca,'diE?t',4);
    irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.95]);
    ylabel(hca,'diE_Z E.B=0')
end
if 1 % angle
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim([diE3t(:,1),angle3t],tint),'k');hold(hca,'on');
    irf_plot(hca,irf_tlim([diE4t(:,1),angle4t],tint),'r');hold(hca,'on'); 
    irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    ylabel('Angle (BtoSpinPlane)')
end

irf_zoom(h,'x',tint)
%% Transform to field aligned coordinate system
t1=[2003 09 29 10 19 27.90];
t2=[2003 09 29 10 19 28.10];
tint=[toepoch(t1) toepoch(t2)];

meanB=mean([irf_tlim(diB1,tint);irf_tlim(diB2,tint);...
        irf_tlim(diB3,tint);irf_tlim(diB4,tint)]);
x=irf_norm(meanB(2:4));    
y=irf_norm(irf_cross(x,irf_cross([1 0 0],x)));
z=irf_norm(irf_cross(x,y));

c_eval('facE?=irf_newxyz(diE?p,x,y,z);');
c_eval('facPos?=irf_newxyz(diPos?,x,y,z);');
%% Plot sc positions
if 0
    irf_units
    figure;
    plot3([mPos1(2)]/Units.RE,[mPos1(3)],[mPos1(4)],'k*'); hold on
    plot3([mPos2(2)]/Units.RE,[mPos2(3)]/Units.RE,[mPos2(4)]/Units.RE,'r*')
    plot3([mPos3(2)]/Units.RE,[mPos3(3)]/Units.RE,[mPos3(4)]/Units.RE,'g*')
    plot3([mPos4(2)]/Units.RE,[mPos4(3)]/Units.RE,[mPos4(4)]/Units.RE,'b*')
end
%% Match
c_eval('mPos?=mean(irf_tlim(facPos?,tint));')
dpar=[mPos1(2)-mPos1(2), mPos1(2)-mPos2(2), mPos1(2)-mPos3(2), mPos1(2)-mPos4(2);...
      mPos2(2)-mPos1(2), mPos2(2)-mPos2(2), mPos2(2)-mPos3(2), mPos2(2)-mPos4(2);...
      mPos3(2)-mPos1(2), mPos3(2)-mPos2(2), mPos3(2)-mPos3(2), mPos3(2)-mPos4(2);...
      mPos4(2)-mPos1(2), mPos4(2)-mPos2(2), mPos4(2)-mPos3(2), mPos4(2)-mPos4(2)];
dperp=[sqrt((mPos1(3)-mPos1(3))^2+(mPos1(4)-mPos1(4))^2),...
       sqrt((mPos1(3)-mPos2(3))^2+(mPos1(4)-mPos2(4))^2),...
       sqrt((mPos1(3)-mPos3(3))^2+(mPos1(4)-mPos3(4))^2),...
       sqrt((mPos1(3)-mPos4(3))^2+(mPos1(4)-mPos4(4))^2);...
       sqrt((mPos2(3)-mPos1(3))^2+(mPos2(4)-mPos1(4))^2),...
       sqrt((mPos2(3)-mPos2(3))^2+(mPos2(4)-mPos2(4))^2),...
       sqrt((mPos2(3)-mPos3(3))^2+(mPos2(4)-mPos3(4))^2),...
       sqrt((mPos2(3)-mPos4(3))^2+(mPos2(4)-mPos4(4))^2);...
       sqrt((mPos3(3)-mPos1(3))^2+(mPos3(4)-mPos1(4))^2),...
       sqrt((mPos3(3)-mPos2(3))^2+(mPos3(4)-mPos2(4))^2),...
       sqrt((mPos3(3)-mPos3(3))^2+(mPos3(4)-mPos3(4))^2),...
       sqrt((mPos3(3)-mPos4(3))^2+(mPos3(4)-mPos4(4))^2);...
       sqrt((mPos4(3)-mPos1(3))^2+(mPos4(4)-mPos1(4))^2),...
       sqrt((mPos4(3)-mPos2(3))^2+(mPos4(4)-mPos2(4))^2),...
       sqrt((mPos4(3)-mPos3(3))^2+(mPos4(4)-mPos3(4))^2),...
       sqrt((mPos4(3)-mPos4(3))^2+(mPos4(4)-mPos4(4))^2)]
       
hh=irf_plot('facE?','comp');

%% Compare to ion moments
sclist=[1 3 4];
c_eval('gseVicod?=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''mat'');',sclist);
c_eval('diVicod?=c_coord_trans(''gse'',''isr2'',gseVicod?,''cl_id'',?);',sclist);
c_eval('facVicod?=irf_newxyz(diVicod?,x,y,z);',sclist);

%% Compare to ion pitch angle distributions
