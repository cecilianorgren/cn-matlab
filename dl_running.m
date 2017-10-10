%% time interval
t1=[2007 08 31 10 18 42.50]; t2=[2007 08 31 10 18 43.00];
t1=[2007 08 31 10 18 42.45]; t2=[2007 08 31 10 18 43.20]; % turb + dl + space 
t1=[2007 08 31 10 18 42.70]; t2=[2007 08 31 10 18 43.00];
%t1=[2007 08 31 10 19 04.50]; t2=[2007 08 31 10 19 05.00]; % lhdw
tint=[toepoch(t1) toepoch(t2)];
sclist=3:4;
clear diE3;
%% load data
if ~exist('diE3','var'); load matlabdiEB; end
%if ~exist('diB3','var'); load matlabB; end

%% zoom in
c_eval('diE?=irf_tlim(diE?,tint);',sclist);
c_eval('diB?=irf_tlim(diB?,tint);',sclist);
diE3=irf_resamp(diE3,diE4);
diB4=irf_resamp(diB4,diE4);
diB3=irf_resamp(diB3,diB4);

%c_eval('gsmE?=irf_tlim(gsmE?,tint);',sclist);
%c_eval('gsmB?=irf_tlim(gsmB?,tint);',sclist);
%% make E either perpendicular field or parallel field
% perp is diE?
c_eval('diE?par=irf_edb(diE?,diB?,70,''Epar'');',sclist);

%% make coordinate system, B, projection of B on spin plane, third
diB0=mean([diB3;diB4]);
z1=[0 0 1]; % z_isr2
x1=irf_norm([diB0(2:3) 0]); % B projection in spin plane
y1=irf_cross(z1,x1); % perpendicular to B in spin plane

%% make field aligned coordinate system, with same y as above
diB0=mean([diB3;diB4]);
z2=irf_norm(diB0(2:4)); % z_B
x2=irf_cross(y1,z2); % perpendicular to B in spin plane
y2=irf_cross(z2,x2); %

%% make field aligned coordinate system, x in propagation direction from tool
% Direction
angles=1:3:360;
f_highpass=7;
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(diB3,diE3,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);
        
diB0=mean([diB3;diB4]);
z3=irf_norm(diB0(2:4)); % z_B
x3=c_coord_trans('gsm','dsi',[tint(1),x(i_dir,:)],'cl_id',3); % this is not exactly perp to z3 since z3 is an average 
                   % between C3 and C4, and x3 is only from one sc, C4 in
                   % this case
x3=x3(2:4);                   
x3=irf_cross(irf_cross(z3,x3),z3); % in supposed propagation direction
y3=irf_cross(z3,x3); % third, current normal direction perhaps

%% coordinate system from article, MVA
gsmx4=-[-0.81 0.54 0.25];
gsmy4=-[-0.09 -0.52 0.85];
gsmz4=[0.59 0.66 0.46];
x4=c_coord_trans('gsm','dsi',[tint(1) gsmx4],'cl_id',3); x4=x4(2:4);
y4=c_coord_trans('gsm','dsi',[tint(1) gsmy4],'cl_id',3); y4=y4(2:4);
z4=c_coord_trans('gsm','dsi',[tint(1) gsmz4],'cl_id',3); z4=z4(2:4);

%% look at coordinate systems
if 1
vec=[x1;y1;z1;x2;y2;z2;x3;y3;z3;x4;y4;z4];
desc={'B proj','perp B in spin plane','z_{isr2}',...
    'x2','perp B in spin plane','B',...
    'v_{tool}','current normal','B',...
    'artx','arty','artz'};
linspec={'b','b','b','r','r','r','g','g','g','k','k','k'};
h=cn_plot3d([0 0 0],vec,desc,linspec);
set(h,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
xlabel('x_{ISR2}'); ylabel('y_{ISR2}'); zlabel('z_{ISR2}')
end

%% turn B and E to new coordinate systems
c_eval('ioB?=irf_newxyz(diB?,x1,y1,z1);',sclist);
c_eval('ioE?=irf_newxyz(diE?,x1,y1,z1);',sclist);
c_eval('ioE?par=irf_newxyz(diE?par,x1,y1,z1);',sclist);
c_eval('facB?=irf_newxyz(diB?,x2,y2,z2);',sclist);
c_eval('facE?=irf_newxyz(diE?,x2,y2,z2);',sclist);
c_eval('facE?par=irf_newxyz(diE?par,x2,y2,z2);',sclist);
c_eval('toB?=irf_newxyz(diB?,x3,y3,z3);',sclist);
c_eval('toE?=irf_newxyz(diE?,x3,y3,z3);',sclist);
c_eval('toE?par=irf_newxyz(diE?par,x3,y3,z3);',sclist);
c_eval('artB?=irf_newxyz(diB?,x4,y4,z4);',sclist);
c_eval('artE?=irf_newxyz(diE?,x4,y4,z4);',sclist);
c_eval('artE?par=irf_newxyz(diE?par,x4,y4,z4);',sclist);

%% plot E and B in coordinate system 1 or 2 / perp or par
if 0
for what=1:8;
switch what
    case 1 % fac par 
        c_eval('E?=facE?par;',sclist);
        c_eval('B?=facB?;',sclist);
        title_str='X=perp to B not in spin plane, Y=perp to B and in spin plane, Z=B \newline reconstructed partly parallel field';
        print_str='fac_par';
    case 2 % fac per
        c_eval('E?=facE?;',sclist);
        c_eval('B?=facB?;',sclist);
        title_str='X=perp to B not in spin plane, Y=perp to B and in spin plane, Z=B \newline reconstructed only perpendicular field';
        print_str='fac_per';
    case 3 % io par
        c_eval('E?=ioE?par;',sclist);
        c_eval('B?=ioB?;',sclist);
        title_str='X=B projection on spin plane, Y=perp to B and in spin plane, Z=spin axis \newline reconstructed partly parallel field';
        print_str='io_par';
    case 4 % io per
        c_eval('E?=ioE?;',sclist);
        c_eval('B?=ioB?;',sclist);
        title_str='X=B projection on spin plane, Y=perp to B and in spin plane, Z=spin axis \newline reconstructed only perpendicular field';
        print_str='io_per';
    case 5 % tool par
        c_eval('E?=toE?par;',sclist);
        c_eval('B?=toB?;',sclist);
        title_str='X=propagation direction from tool, Y=current normal?, Z=B \newline reconstructed partly parallel field';
        print_str='to_par';
    case 6 % tool per
        c_eval('E?=toE?;',sclist);
        c_eval('B?=toB?;',sclist);
        title_str='X=propagation direction from tool, Y=current normal?, Z=B \newline reconstructed only perpendicular field';
        print_str='to_per';
    case 7 % isr2 par
        c_eval('E?=diE?par;',sclist);
        c_eval('B?=diB?;',sclist);
        title_str='ISR2 \newline reconstructed partly parallel field';
        print_str='di_par';
    case 8 % isr2 per
        c_eval('E?=diE?;',sclist);
        c_eval('B?=diB?;',sclist);
        title_str='ISR2 \newline reconstructed only perpendicular field';
        print_str='di_per';
end       

fig=figure(99);
set(fig,'position',[560 44 600 900])
setupfigure;
h=irf_plot(6);
isub=1;
if 1 % B_X
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{B3(:,[1 2]),B4(:,[1 2])},'comp')
    ylabel(hca,'B_X')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Y
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{B3(:,[1 3]),B4(:,[1 3])},'comp')
    ylabel(hca,'B_Y')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Z
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{B3(:,[1 4]),B4(:,[1 4])},'comp')
    ylabel(hca,'B_Z')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % E_X
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{E3(:,[1 2]),E4(:,[1 2])},'comp')
    ylabel(hca,'E_X')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % E_Y
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{E3(:,[1 3]),E4(:,[1 3])},'comp')
    ylabel(hca,'E_Y')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % E_Z
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{E3(:,[1 4]),E4(:,[1 4])},'comp')
    ylabel('E_Z')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
irf_zoom(h,'x',tint)
title(h(1),title_str)
eval(['print -dpng /Users/Cecilia/LHDW/Pics/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_',print_str,'.gif']);
end
end

%% look at the topology of the E and B fields
% use cn_plot3d and look to dl_b

% load positions of satellites
c_eval('gsePos?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);         
c_eval('gsmPos?=c_coord_trans(''GSE'',''GSM'',gsePos?,''cl_id'',?);',sclist);
%c_eval('diPos?=c_coord_trans(''GSE'',''ISR2'',gsePos?,''cl_id'',4);',sclist);
clear gsePos3 gsePos4;

% zoom in on instantaneous position and average 
c_eval('gsmPos?=irf_tlim(gsmPos?,tint);',sclist)
c_eval('gsmPos?av=mean(gsmPos?(:,2:end));',sclist) % km 
%c_eval('diPos?=irf_tlim(diPos?,tint);',sclist)
%c_eval('diPos?av=mean(diPos?(:,2:end));',sclist) % km 

% get zoomed in relative positions
%R0=(diPos3av+diPos4av)/2;
%c_eval('diPos?avrel=diPos?av-R0;',sclist);
%c_eval('R?=diPos?avrel;',sclist);

gsmR0=(gsmPos3av+gsmPos4av)/2;
c_eval('gsmPos?avrel=gsmPos?av-gsmR0;',sclist);
c_eval('gsmR?=gsmPos?avrel;',sclist);
c_eval('R?=c_coord_trans(''gsm'',''dsi'',[tint(1) gsmR?],''cl_id'',3);',[0 sclist]);
c_eval('R?=R?(2:4);',sclist);

% get positions in different coordinate systems
c_eval('ioR?=irf_newxyz(R?,x1,y1,z1);',[0 sclist]);
c_eval('facR?=irf_newxyz(R?,x2,y2,z2);',[0 sclist]);
c_eval('toR?=irf_newxyz(R?,x3,y3,z3);',[0 sclist]);
c_eval('artR?=irf_newxyz(R?,x4,y4,z4);',[0 sclist]);

% velocity and diff distance vector in relevant coordinate system
v=860;
vamp=[v 0 0]; % x y z km/s in any coordinate system where x is propagation direction
dt=diff(diE3(1:2,1));
dr=vamp*dt; % km
dt_match=0.0376; % ms
%dx_match4=dt_match*to_vamp;
dx_match3=[0 0 0];
dx_match4=dx_match3;
dt_laps=dt_match/dt;
c_eval('');



% get subsequent positions of C3 C4
RO_vec=bsxfun(@times,bsxfun(@minus,diE3(:,1),diE3(1,1)),-vamp);
c_eval('O?=bsxfun(@plus,RO_vec,toR?-dx_match?);',sclist);

% calculate parallel current density j_z=mu0*(dBy/dx-dBx/dy)
units=irf_units;
toj_z=[toB3(:,1) units.mu0*((toB4(:,3)-toB3(:,3))/(toR4(1)-toR3(1))-(toB4(:,2)-toB3(:,2))/(toR4(2)-toR3(2)))];
facj_z=[facB3(:,1) units.mu0*((facB4(:,3)-facB3(:,3))/(facR4(1)-facR3(1))-(facB4(:,2)-facB3(:,2))/(facR4(2)-facR3(2)))];
artj_z=[artB3(:,1) units.mu0*((artB4(:,3)-artB3(:,3))/(artR4(1)-artR3(1))-(artB4(:,2)-artB3(:,2))/(artR4(2)-artR3(2)))];

% filter B and E
f_filt=0.1;
c_eval('artB?ac=irf_filt(artB?,f_filt,0,[],5);',sclist);
c_eval('artE?ac=irf_filt(artE?,f_filt,0,[],5);',sclist);

c_eval('toB?ac=irf_filt(toB?,f_filt,0,[],5);',sclist);
c_eval('toE?ac=irf_filt(toE?,f_filt,0,[],5);',sclist);
c_eval('toE?ac=toE?;',sclist);

%% now we plot magnetic field with Epar
%text(toR4(1),toR4(2),toR4(3),'C4'); hold on;
figure
quiver3(O3(:,1),O3(:,2),O3(:,3),toB3ac(:,2),toB3ac(:,3),toB3ac(:,4),'r'); hold on;
quiver3(O4(:,1),O4(:,2),O4(:,3),toB4ac(:,2),toB4ac(:,3),toB4ac(:,4),'r');
quiver3(O3(:,1),O3(:,2),O3(:,3),toE3(:,2),toE3(:,3),toE3(:,4),'b');
quiver3(O4(:,1),O4(:,2),O4(:,3),toE4(:,2),toE4(:,3),toE4(:,4),'b');
quiver3(O3(:,1),O3(:,2),O3(:,3),toE3par(:,2),toE3par(:,3),toE3par(:,4),'g');
quiver3(O4(:,1),O4(:,2),O4(:,3),toE4par(:,2),toE4par(:,3),toE4par(:,4),'g');
c_eval('text(toR?(1),toR?(2),toR?(3),''C?''); hold on;',sclist);
axis equal
xlabel('x');ylabel('y');zlabel('z')

%% now we plot
quiver3(O3(:,1),O3(:,2),O3(:,3),toE3par(:,2),toE3par(:,3),toE3par(:,4)); hold on;
quiver3(O4(:,1),O4(:,2),O4(:,3),toE4par(:,2),toE4par(:,3),toE4par(:,4));
quiver3(O3(:,1),O3(:,2),O3(:,3),toE3(:,2),toE3(:,3),toE3(:,4));
quiver3(O4(:,1),O4(:,2),O4(:,3),toE4(:,2),toE4(:,3),toE4(:,4));
axis equal
xlabel('k_tool'); ylabel('third'); zlabel('B')

%% now we plot magnetic field with Eperp
% first must filter B
c_eval('toB?ac=irf_filt(toB?,1,0,[],5);',sclist);
quiver3(O3(:,1),O3(:,2),O3(:,3),toB3ac(:,2),toB3ac(:,3),toB3ac(:,4),'r'); hold on;
quiver3(O4(:,1),O4(:,2),O4(:,3),toB4ac(:,2),toB4ac(:,3),toB4ac(:,4),'r');
quiver3(O3(:,1),O3(:,2),O3(:,3),toE3(:,2),toE3(:,3),toE3(:,4),'b');
quiver3(O4(:,1),O4(:,2),O4(:,3),toE4(:,2),toE4(:,3),toE4(:,4),'b');
axis equal

%% make gif with parallel propagation and a plane moving along to see projection of fields.
% use dt_laps to cut away the first times on C4, and last on C3.
E4cut=toE4par(round(dt_laps+2):end,:);
E3cut=toE3par(1:(end-round(dt_laps-1)),:);
B4cut=irf_norm(toB4ac(round(dt_laps+2):end,:));
B3cut=irf_norm(toB3ac(1:(end-round(dt_laps-1)),:));

E4cut_norm=irf_norm(toE4(round(dt_laps+2):end,:));
E3cut_norm=irf_norm(toE3(1:(end-round(dt_laps-1)),:));

%c_eval('E?cut_norm=irf_norm(E?cut);',sclist)

B4cutdc=irf_norm(toB4(round(dt_laps+2):end,:));
B3cutdc=irf_norm(toB3(1:(end-round(dt_laps-1)),:));
%plot(E3cut(:,4));hold on;
%plot(E4cut(:,4))
%%


fig=figure(35);
set(fig,'position',  [933      109     450       925 ]);
h(1)=axes('position',[0.070    0.64    0.8750    0.34]);
h(2)=axes('position',[0.070    0.51    0.8750    0.09]);
h(3)=axes('position',[0.070    0.42    0.8750    0.09]);
h(4)=axes('position',[0.070    0.33    0.8750    0.09]);
h(5)=axes('position',[0.070    0.24    0.8750    0.09]);
h(6)=axes('position',[0.070    0.15    0.8750    0.09]);
h(7)=axes('position',[0.070    0.06    0.8750    0.09]);
%
setupfigure
cC3=[0 1 0.3];
cC4=[0 0.3 1];
scale=10;
hold3='on'; hold4='off';
c_eval('text(toR?(1),toR?(2),toR?(3),''C?''); hold on;',sclist);
axis(h(1),'equal')
for k=1:size(E3cut,1)-2
    ind=k;    
    %c_eval('quiver(h(1),toR?(1),toR?(2),B?cut(k-2,2)*scale,B?cut(k-2,3)*scale,''color'',[0.5 0.5 1]); hold(h(1),hold3);',sclist);
    %c_eval('quiver(h(1),toR?(1),toR?(2),B?cut(k-1,2)*scale,B?cut(k-1,3)*scale,''color'',[0.5 0.5 1]); hold(h(1),hold3);',sclist);
    c_eval('quiver(h(1),toR?(1),toR?(2),E?cut_norm(k,2)*scale,E?cut_norm(k,3)*scale,''color'',[1 0 1]); hold(h(1),hold3);',sclist);
    c_eval('quiver(h(1),toR?(1),toR?(2),B?cut(k,2)*scale,B?cut(k,3)*scale,''color'',[0 0 1]); hold(h(1),hold3);',sclist);
    c_eval('plot(h(1),toR?(1)+B?cut(1:k,2)*scale,toR?(2)+B?cut(1:k,3)*scale,''color'',cC?)',sclist);
    axis(h(1),'equal');hold(h(1),'off')
    %quiver(h(1),toR4(1),toR4(2),B4cut(k,2)*scale,B4cut(k,3)*scale); hold(h(1),'off');    
    set(h(1),'xlim',(abs(toR3(1))+scale)*[-1 1],'ylim',(abs(toR3(2))+scale)*[-1 1])
    
    %pause(0.2);
    set(h(2:7),'ColorOrder',[cC3;cC4])
    c_eval('irf_plot(h(?),{E3cut(:,[1 ?]),E4cut(:,[1 ?])},''comp''); hold(h(?),''on'');',2:4);
    c_eval('irf_plot(h(?),E3cut(k,[1 ?]),''ro'');irf_plot(h(?),E4cut(k,[1 ?]),''ro'');',2:4);
    c_eval('irf_plot(h(?+3),{B3cutdc(:,[1 ?]),B4cutdc(:,[1 ?])},''comp''); hold(h(?+3),''on'');',2:4);
    c_eval('irf_plot(h(?+3),B3cutdc(k,[1 ?]),''ro'');irf_plot(h(?+3),B4cutdc(k,[1 ?]),''ro'');',2:4);
    c_eval('hold(h(?),''off'')',2:7);
    c_eval('grid(h(?),''off'')',2:7);
    irf_zoom(h(2:7),'x',[min([E3cut(1,1) E4cut(1,1)]) max([E3cut(end,1) E4cut(end,1)])]);
     
    if 1
% collect frames    
    f=getframe(fig);
    A(:,ind)=f;
    if k==1, % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,size(E3cut,1)) = 0;
    else
        im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
    end
    end
end
imwrite(im,map,['/Users/Cecilia/DL/Pics/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_C34_Bperp.gif'],'DelayTime',0.01,'LoopCount',inf);

%%
%desc={'B proj','perp B in spin plane','z_{isr2}',...
%    'x2','perp B in spin plane','B',...
%    'v_{tool}','current normal','B'};
desc={' '};
linspec={'g','b'};
plot(1);
h=gca;
cn_plot3d(h,O3,diE3(:,2:4),' ','b');
%cn_plot3(diE3(:,2:4),diE4(:,2:4),)
set(h,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
xlabel('x_{ISR2}'); ylabel('y_{ISR2}'); zlabel('z_{ISR2}')


%% load/construct polar/cross-section particle data
%c_eval('distr?_3dxph=cn_construct_distribution_data_combined(''C?_CP_PEA_3DXPH_PSD'');',sclist);
%c_eval('distr?_pitch_3dxh=cn_construct_distribution_data_combined(''C?_CP_PEA_PITCH_3DXH_PSD'');',sclist);
%c_eval('distr?_pitch_3drl=cn_construct_distribution_data_combined(''C3_CP_PEA_PITCH_3DRL_PSD'');',3);
%c_eval('distr?_pitch_full_leea=cn_construct_distribution_data_combined(''C3_CP_PEA_PITCH_FULL_PSD'',''LEEA'');',3);
%c_eval('distr?_pitch_full_heea=cn_construct_distribution_data_combined(''C3_CP_PEA_PITCH_FULL_PSD'',''HEEA'');',3);
load distributions

%% plot particles
fig=figure(54);
set(fig,'position',[1 -21 1676 955]);
dx=0.18; h=0.35; w=0.15;
for k=1:5; h(k)=axes('position',[0.05+(k-1)*dx 0.56 0.14 0.35]); end
for k=6:10; h(k)=axes('position',[0.05+(k-5-1)*dx 0.10 0.14 0.35]); end

dt=0.13;
for k=1:5
    tint_k=tint(1)-0.2+(k-1)*dt;
    cn_plot_distribution_function_combined(h(k),tint_k,'cross-section',distr3_pitch_full_leea);
    hold(h(k),'on')
    cn_plot_distribution_function_combined(h(k+5),tint_k,'cross-section',distr3_pitch_full_heea);
    hold(h(k),'off');
end

set(fig,'PaperPositionMode','auto');
%% load/construct polar/cross-section particle data
c_eval('psd?pol=cn_construct_distribution_data(''C?_CP_PEA_3DXPH_PSD'',''polar'');',sclist);
c_eval('psd?cs=cn_construct_distribution_data(''C?_CP_PEA_3DXPH_PSD'',''cross-section'');',sclist);
c_eval('psd?pol_pf=cn_construct_distribution_data(''C?_CP_PEA_PITCH_FULL_PSD'',''polar'');',sclist);
%% many times in different subplots, c3
for k=1:10; h(k)=subplot(2,5,k); end
isub=1;
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3cs); 
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-5-1)*0.125,psd3pol_pf);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3pol_pf);
end
if 1
    hca=h(isub); isub=isub+1
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3pol_pf);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3pol_pf);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-1)*0.125,psd3pol_pf);
end
%% many times in different subplots, c3
fig=figure(54);
set(fig,'position',[1         -21        1676         955])
for k=1:10; h(k)=subplot(2,5,k); end
isub=1;
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-2)*0.125,psd4cs); 
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-2)*0.125,psd4cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-2)*0.125,psd4cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-2)*0.125,psd4cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-2)*0.125,psd4cs);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-5-2)*0.125,psd4pol);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-5-2)*0.125,psd4pol);
end
if 1
    hca=h(isub); isub=isub+1
    cn_plot_distribution_function(hca,tint(1)+(isub-5-2)*0.125,psd4pol);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-5-2)*0.125,psd4pol);
end
if 1
    hca=h(isub); isub=isub+1;
    cn_plot_distribution_function(hca,tint(1)+(isub-5-2)*0.125,psd4pol);
end
%% compare new and old distribution_dunction... cross-section
h1=subplot(1,2,1);
c_caa_distribution_function(h1,'tint',[tint(1) tint(2)],'C4_CP_PEA_3DXPH_PSD','cross-section');
h2=subplot(1,2,2);
cn_plot_distribution_function(h2,tint,psd4cs);
setupfigure
%% ok
%distr_3dxh_3=cn_construct_distribution_data_combined('C3_CP_PEA_PITCH_3DXH_PSD');
distr_3dxph_4=cn_construct_distribution_data_combined('C4_CP_PEA_3DXPH_PSD');
for k=1:6
    hca=subplot(2,3,k);
    cn_plot_distribution_function_combined(hca,tint(1)-0.5+(k-1)*0.130,'polar',distr_3dxh_3);
end
%% compare new and old distribution_dunction... polar
%c_eval('psd?pol=cn_construct_distribution_data_combined(''C?_CP_PEA_3DXPH_PSD'');',sclist);
h1=subplot(1,3,1);
%cn_plot_distribution_function(h1,[tint(2) tint(2)+0.13],psd4pol);
cn_plot_distribution_function_combined(h1,[tint(1) tint(1)+0.13],'cross-section',distr4_3dxph);
h2=subplot(1,3,2);
cn_plot_distribution_function_combined(h2,[tint(1) tint(1)+0.13],'polar',distr_3dxph_4);
h3=subplot(1,3,3);
cn_plot_distribution_function_combined(h3,[tint(1) tint(1)+0.13],'cross-section',distr_3dxph_4);
%c_caa_distribution_function(h3,'tint',[tint(1) tint(2)],'C4_CP_PEA_3DXPH_PSD','polar');
setupfigure