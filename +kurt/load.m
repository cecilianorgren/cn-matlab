%% Load data to investigate Kurt's question
cd /Users/Cecilia/Data/BM/20070831

% load magnetic field 
c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
load mBS; tint=[dBS3(1,1) dBS3(end,1)];
c_eval('gseB?staff=c_coord_trans(''DSC'',''gse'',dBS?,''cl_id'',?);',3:4);
c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?staff);',3:4)
offset=cn_offset(3,gseB3fgm,4,gseB4fgm,'gsm');  % whole BM interval
offset2=cn_offset(...
    3,irf_tlim(gseB3fgm,toepoch([2007 08 31 10 10 0]),toepoch([2007 08 31 10 13 0])),...
    4,irf_tlim(gseB4fgm,toepoch([2007 08 31 10 10 0]),toepoch([2007 08 31 10 13 0])),...
    'gsm');                                     % out in lobes
gseB3fgm_offs=cn_apply_offset(3,gseB3fgm,'gsm',offset2);
gseB3_offs=c_fgm_staff_combine(gseB3fgm_offs(:,1:4),gseB3staff);

gsedB = irf_add(-1,gseB3,1,gseB4); % B4-B3
gsedB_offs = irf_add(-1,gseB3_offs,1,gseB4); % B4-B3offs

% load electric field
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);
c_eval('gseE?=c_coord_trans(''dsi'',''gse'',diE?,''CL_ID'',?);',3:4);

%% Create field aligned coordinate system
tint = toepoch([2007 08 31 10 18 57;2007 08 31 10 19 13])';
n_hat = -[0.0856    0.9293   -0.3594]; % gse
n_hat = -[0.2315    0.8372   -0.4955]; % cn_tool
% n_hat = [-0.09 -0.89 0.44]; % from article, gsm -> gse
b_hat = irf_norm(irf_add(0.5,irf_tlim(gseB3,tint),0.5,irf_tlim(gseB4,tint))); % gse
n_hat = irf_norm(irf_resamp([tint(1) n_hat],b_hat));
v_hat = irf_norm(irf_cross(b_hat,n_hat));

% x: n_hat
% y: v_hat
% z: b_hat

c_eval('b1? = irf_dot(n_hat,gseB?);',3:4) % normal to boundary layer
c_eval('b2? = irf_dot(v_hat,gseB?);',3:4) % propagation direction
c_eval('b3? = irf_dot(b_hat,gseB?);',3:4) % magnetic field direction
c_eval('facB? = [b1? b2?(:,2) b3?(:,2)];',3:4); 
clear b13 b23 b33 b14 b24 b34
c_eval('e1? = irf_dot(n_hat,gseE?);',3:4) % normal to boundary layer
c_eval('e2? = irf_dot(v_hat,gseE?);',3:4) % propagation direction
c_eval('e3? = irf_dot(b_hat,gseE?);',3:4) % magnetic field direction
c_eval('facE? = [e1? e2?(:,2) e3?(:,2)];',3:4); 
clear e13 e23 e33 e14 e24 e34
db1 = irf_dot(n_hat,gsedB); % normal to boundary layer
db2 = irf_dot(v_hat,gsedB); % propagation direction
db3 = irf_dot(b_hat,gsedB); % magnetic field direction
facdB = [db1 db2(:,2) db3(:,2)];
clear db13 db23 db33 db14 db24 db34
db1 = irf_dot(n_hat,gsedB_offs); % normal to boundary layer
db2 = irf_dot(v_hat,gsedB_offs); % propagation direction
db3 = irf_dot(b_hat,gsedB_offs); % magnetic field direction
facdB_offs = [db1 db2(:,2) db3(:,2)];
clear db13 db23 db33 db14 db24 db34
c_eval('sb1? = irf_dot(n_hat,gseB?staff);',3:4) % normal to boundary layer
c_eval('sb2? = irf_dot(v_hat,gseB?staff);',3:4) % propagation direction
c_eval('sb3? = irf_dot(b_hat,gseB?staff);',3:4) % magnetic field direction
c_eval('facB?staff = [sb1? sb2?(:,2) sb3?(:,2)];',3:4); 
clear sb13 sb23 sb33 sb14 sb24 sb34

%% Integrate electric field
tref = toepoch([2007 08 31 10 19 4.8]);
v = 1000; % km/s 
f_filt = 10;
c_eval('facEfilt? = irf_filt(facE?,5,0,450,3);',3:4)
c_eval('P? = irf_integrate(irf_tlim(facE?(:,[1 3]),tint),tref);',3:4);
c_eval('P? = irf_tappl(P?,''*1000'');',3:4)
c_eval('P? = irf_filt(P?,f_filt,0,450,3);',3:4)
c_eval('facBfilt? = irf_filt(facB?,f_filt,0,450,3);',3:4)
%% Plot magnetic field and electric field
if 1
h=irf_plot(7); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(facE3,tint)); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'E_{FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,facB3); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,facdB); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B4-B3 [FAC]') 
hca=h(isub);isub=isub+1;irf_plot(hca,facdB_offs); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B4-B3_{offs} [FAC]') 
hca=h(isub);isub=isub+1;irf_plot(hca,{P3,P4},'comp'); irf_legend(hca,{'C3','C4'},[0.95 0.92]); ylabel(hca,'\phi'); irf_zoom(hca,'y',550*[-1 1])
hca=h(isub);isub=isub+1;irf_plot(hca,facBfilt3(:,[1 4])); ylabel(hca,'B_{||,FGM+STAFF,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,facBfilt3); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
end

%% Plot perdpendicular field to see structure
% load satellite positions
if ~exist('facR3','var')
    c_eval('gseR?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
    c_eval('r1? = irf_dot(n_hat,gseR?);',3:4) % normal to boundary layer
    c_eval('r2? = irf_dot(v_hat,gseR?);',3:4) % propagation direction
    c_eval('r3? = irf_dot(b_hat,gseR?);',3:4) % magnetic field direction
    c_eval('facR? = [r1? r2?(:,2) r3?(:,2)];',3:4); 
    clear r13 r23 r33 r14 r24 r34
end

propdir = 'rightleft';
tint = toepoch([2007 08 31 10 19 04.5; 2007 08 31 10 19 05.2])';
%tint = toepoch([2007 08 31 10 19 05.5 ;2007 08 31 10 19 05.9])';
v = -800;
scaleE = 12*2;
scaleB = 3*2;
kurt.ebplot(irf_tlim(facEfilt3,tint),...
            irf_tlim(facEfilt4,tint),...
            irf_tlim(facBfilt3,tint),...
            irf_tlim(facBfilt4,tint),...
            irf_tlim(facR3,tint),...
            irf_tlim(facR4,tint),...
            v,propdir,scaleE,scaleB);
        
%% Check the angle between dB and dE 
c_eval('normB? = irf_norm(facBfilt?);',3:4)
normB3(:,3) = -normB3(:,3);
c_eval('dot? = irf_tlim(irf_dot(irf_norm(facEfilt?) ,normB?),tint);',3:4)
c_eval('dot? = irf_tlim(irf_dot(irf_norm(facEfilt?),irf_norm(facBfilt?)),tint);',3:4)
c_eval('deg? = [dot?(:,1) acosd(dot?(:,2))];',3:4)

if 1 % Plot current and etc
h=irf_plot(5); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,{facEfilt3(:,[1 2]),facEfilt4(:,[1 2])},'comp'); ylabel(hca,'E_n [mV/m]'); irf_legend(hca,{'C3','C4'},[0.02 0.9])
hca=h(isub);isub=isub+1;irf_plot(hca,{facEfilt3(:,[1 3]),facEfilt4(:,[1 3])},'comp'); ylabel(hca,'E_k [mV/m]'); irf_legend(hca,{'C3','C4'},[0.02 0.9])
hca=h(isub);isub=isub+1;irf_plot(hca,{facBfilt3(:,[1 2]),facBfilt4(:,[1 2])},'comp'); ylabel(hca,'B_n [nT]'); irf_legend(hca,{'C3','C4'},[0.02 0.9])
hca=h(isub);isub=isub+1;irf_plot(hca,{irf_tappl(facBfilt3(:,[1 3]),'*-1'),facBfilt4(:,[1 3])},'comp'); ylabel(hca,'B_k [nT]'); irf_legend(hca,{'C3','C4'},[0.02 0.9])
hca=h(isub);isub=isub+1;irf_plot(hca,{deg3,deg4},'comp'); ylabel(hca,'Angle between B and E'); irf_legend(hca,{'C3','C4'},[0.02 0.9])
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
end
%% Calculate parallel current from knowing dB along two lines
% use c_4_j.m:
%c_4_j  Calculate current from using 4 spacecraft technique
%   in addition one can obtain average magnetic field and jxB values
% 
%   [j,divB,B,jxB,curvature,divTshear,divPb] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4)
%   [j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?')
%   Estimate also divergence B as the error estimate

tint = toepoch([2007 08 31 10 19 04.5 ;2007 08 31 10 19 05.2])';
% time transformed into length 
v = -1100; % km/s

% try first by making a square box.
% the fixed direction is the normal direction, x.
meanR3 = cn.mean(irf_tlim(facR3,tint),1);
meanR4 = cn.mean(irf_tlim(facR4,tint),1);
Rn = meanR4(1) - meanR3(1);
Rk = meanR4(2) - meanR3(2);

tdiff = Rn/v;
tdiff2 = Rk/v;

facR3a = facR3; facR3a(:,3) = facR3a(:,3) + 0*Rk;
facR4a = facR4; facR4a(:,3) = facR4a(:,3) + 0*Rk;
facR3b = facR3; facR3b(:,3) = facR3b(:,3) + Rn; %facR3b(:,1) = facR3b(:,1) + tdiff;
facR4b = facR4; facR4b(:,3) = facR4b(:,3) + Rn; %facR4b(:,1) = facR4b(:,1) + tdiff;

% adjust time series and make new time shifted magnetic field vectors
facB3b = facBfilt3; facB3b(:,1) = facB3b(:,1) + tdiff;
facB4b = facBfilt4; facB4b(:,1) = facB4b(:,1) + tdiff;
facB3a = facBfilt3; facB3a(:,1) = facB3a(:,1) + 0*tdiff2;
facB4a = facBfilt4; facB4a(:,1) = facB4a(:,1) + 0*tdiff2;

[j jrot] = kurt.j(facR3b,facR4b,facR3a,facR4a,facB3b,facB4b,facB3a,facB4a,155);
%[j,divB,B,jxB,divTshear,divPb] = c_4_j(facR3b,facR4b,facR3a,facR4a,facB3b,facB4b,facB3a,facB4a);

if 0 % Plot current and etc
h=irf_plot(7); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(j,tint)); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'j_{FAC} [A?]'); irf_zoom(hca,'y',[-1 1]*1e-5)
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(j(:,[1 4]),tint)); ylabel(hca,'j_{||,FAC} [A?]'); irf_zoom(hca,'y',[-1 1]*1e-5)
hca=h(isub);isub=isub+1;irf_plot(hca,facB3); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
%hca=h(isub);isub=isub+1;irf_plot(hca,facdB); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B4-B3 [FAC]') 
hca=h(isub);isub=isub+1;irf_plot(hca,facdB_offs); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B4-B3_{offs} [FAC]') 
%hca=h(isub);isub=isub+1;irf_plot(hca,{P3,P4},'comp'); irf_legend(hca,{'C3','C4'},[0.95 0.92]); ylabel(hca,'\phi'); irf_zoom(hca,'y',550*[-1 1])
hca=h(isub);isub=isub+1;irf_plot(hca,facBfilt3(:,[1 4])); ylabel(hca,'B_{||,FGM+STAFF,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,facBfilt3); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facB3b,facB4b,facB3,facB4},'comp'); irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facR3b,facR4b,facR3,facR4},'comp'); irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.92]); ylabel(hca,'R_{FAC}')
irf_zoom(h,'x',tint)
end

if 1 % Plot current and etc
h=irf_plot(4); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(j,tint)); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'j [A?]'); %irf_zoom(hca,'y',[-1 1]*1e-5)
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(jrot,tint)); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'j_{rot} [A?]'); %irf_zoom(hca,'y',[-1 1]*1e-5)
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(irf_tappl(j(:,[1 4]),'*1e9'),tint)); ylabel(hca,'j_{||,FAC} [nA?]'); %irf_zoom(hca,'y',[-1 1]*1e-5)
%hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(divB,tint)); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'div B'); irf_zoom(hca,'y',[-1 1]*1e-5)
%hca=h(isub);isub=isub+1; irf_plot(hca,facdB_offs); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B4-B3_{offs}') 
hca=h(isub);isub=isub+1;irf_plot(hca,facBfilt3); irf_legend(hca,{'n','v','b'},[0.95 0.92]); ylabel(hca,'B')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facB3b(:,[1 2]),facB4b(:,[1 2]),facB3(:,[1 2]),facB4(:,[1 2])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'B_{n}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facB3b(:,[1 3]),facB4b(:,[1 3]),facB3(:,[1 3]),facB4(:,[1 3])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'B_{k}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facB3b(:,[1 4]),facB4b(:,[1 4]),facB3(:,[1 4]),facB4(:,[1 4])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'B_{b}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facR3b(:,[1 2]),facR4b(:,[1 2]),facR3(:,[1 2]),facR4(:,[1 2])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'R_{n}'); irf_zoom(hca,'y');
%hca=h(isub);isub=isub+1;irf_plot(hca,{facR3b(:,[1 3]),facR4b(:,[1 3]),facR3(:,[1 3]),facR4(:,[1 3])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'R_{k}'); irf_zoom(hca,'y');
%hca=h(isub);isub=isub+1;irf_plot(hca,{facR3b(:,[1 4]),facR4b(:,[1 4]),facR3(:,[1 4]),facR4(:,[1 4])},'comp'); irf_legend(hca,{'C3b','C4b','C3','C4'},[0.95 0.92]); ylabel(hca,'R_{b}'); irf_zoom(hca,'y');
title(h(1),'FAC')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facB3b,facB4b,facB3,facB4},'comp'); irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.92]); ylabel(hca,'B_{FGM+STAFF,FAC}')
%hca=h(isub);isub=isub+1;irf_plot(hca,{facR3b,facR4b,facR3,facR4},'comp'); irf_legend(hca,{'C1','C2','C3','C4'},[0.95 0.92]); ylabel(hca,'R_{FAC}')
irf_zoom(h,'x',tint)
end
%%







