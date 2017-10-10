% investigate how waveB is round dl.
t1=[2007 08 31 10 18 42.50]; t2=[2007 08 31 10 18 43.00];
tint=[toepoch(t1) toepoch(t2)];
sclist=3:4;

% load data
if ~exist('diE3','var'); load matlabE; end
if ~exist('diB3','var'); load matlabB; end
%c_eval('diB?=c_coord_trans(''GSM'',''ISR2'',gsmB?,''cl_id'',?);',sclist);

% zoom in
c_eval('E?=irf_tlim(diE?,tint);',sclist);
c_eval('B?=irf_tlim(diB?,tint);',sclist);

%c_eval('gsm_B?=irf_tlim(gsmB?,tint);',sclist);
%gsmB0=mean([gsm_B3(:,2:4); gsm_B4(:,2:4)],1);
B0=mean([B3(:,2:4); B4(:,2:4)],1);
c_eval('B0?=mean(B?(:,2:4),1);',sclist);

% coordinate system details
normal=-[0.0856 -0.52 0.85]; % current sheet normal in GSM
di_normal=c_coord_trans('GSM','ISR2',[B3(1,1) normal],'cl_id',3); % current sheet normal in ISR2
di_normal=di_normal(2:4);
z=irf_norm(B0); % B in ISR2
y=irf_norm(irf_cross(irf_cross(z,di_normal),z)); % current sheet normal in ISR2
x=irf_cross(y,z); % current direction in ISR2
inspinplane=irf_cross(z,[0 0 1]); % measured component in spin plane, ISR2


c_eval('z?=mean(B?(:,2:4),1);',sclist);


% construct field aligned coordinate system
%facB3=irf_newxyz(B3,x,y,z);
%facB4=irf_newxyz(B4,x,y,z);
c_eval('facB?=irf_newxyz(B?,x,y,z);',sclist);
c_eval('facB??=irf_newxyz(B?,x,y,z?);',sclist);
%%
figure(69);
setupfigure;
h=irf_plot({facB3,facB4},'comp');
ylabel(h(1),'B_X (current)')
ylabel(h(2),'B_Y (current normal)')
ylabel(h(3),'B_Z (B_0)')
ylabel(h(4),'|B|')
for k=1:4; 
    irf_legend(h(k),{'C3','C4'},[0.95 0.95]);  
    irf_zoom(h,'y')
end

%% make correlation and find out angle between the B field at the two satellites
%[t02DC corrxDC corr_maxDC]=cn_xcorr(detrend(facB4(:,2:4),'linear'),detrend(facB3(:,2:4),'linear'),window,'x');
B4=[facB4(:,1) detrend(facB4(:,2:4),'linear')];
B3=[facB3(:,1) detrend(facB3(:,2:4),'linear')];
[corr lags]=xcorr(detrend(facB4(:,2),'linear'),detrend(facB3(:,2),'linear'));
dlags=lags(find(corr==max(corr)));
dt=dlags/450;
dot34=irf_dot([B3(dlags+1:end,:)],B4);
dot34_per=irf_dot([B3(dlags+1:end,1:3) B3(dlags+1:end,1)*0],[B4(:,1:3) B4(:,1)*1]);
angle=[dot34(:,1) acosd(dot34(:,2))];
angle_per=[dot34_per(:,1) acosd(dot34_per(:,2))];
%% draw topology, detrend field
cn_plot3([facB3(:,1) detrend(facB3(:,2:4))],[facB4(:,1) detrend(facB4(:,2:4))])
%%  draw topology, filter field
cn_plot3(irf_filt(facB3,1,180,450,5),irf_filt(facB4,1,180,450,5))

