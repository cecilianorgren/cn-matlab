t1=[2007 08 31 10 19 04.50];  t2=[2007 08 31 10 19 05.10]; 
n_hat=-[0.0856    0.9293   -0.3594];
%n_hat=-[0.2315    0.8372   -0.4955]; cn_tool
t1=[2007 08 31 10 19 07.20];  t2=[2007 08 31 10 19 07.50];
t1=[2007 08 31 10 19 05.50];  t2=[2007 08 31 10 19 05.90];

n_hat=-[0.0856    0.9293   -0.3594];
t1    =[2007 08 31 10 19 05.50];
t2    =[2007 08 31 10 19 05.90];
% DL
n_hat=-[0.0856    0.9293   -0.3594];
t1    =[2007 08 31 10 18 42.50];
t2    =[2007 08 31 10 18 43.00];
%%
n_hat=-[ 0.113  0.913 -0.391]; % gse
%%
t1str=[num2str(t1(1)),'0',num2str(t1(2)),'',num2str(t1(3)),...
            '',num2str(t1(4)),'',num2str(t1(5)),'',num2str(t1(6))];
t2str=[num2str(t2(1)),'0',num2str(t2(2)),'',num2str(t2(3)),...
            '',num2str(t2(4)),'',num2str(t2(5)),'',num2str(t2(6))];

% Taking average density during time interval
c_eval('peaNe?z=cn_toepoch(t1,t2,peaNe?);',3:4);
%c_eval('peaNe?av=sum(peaNe?z(:,2),1)/size(peaNe?z,1);',3:4);
c_eval('Ne?z=irf_resamp(peaNe?z,[peaNe?(1,1):1/450:peaNe?z(end,1)]'');',3:4)
c_eval('Ne?z=cn_toepoch(t1,t2,Ne?z);',3:4);
Neav=mean([Ne3z(:,2); Ne4z(:,2)]);
%peaNeav=(Ne3av+Ne4av)/2;
%Neav=peaNeav;

Nizoom=cn_toepoch(t1,t2,codNi4);
c_eval('Ni?z=irf_resamp(Nizoom,[Nizoom(1,1):1/450:Nizoom(end,1)]'');',4)
Ni4z=cn_toepoch(t1,t2,Ni4z);
Niav=mean(Ni4z(:,2));

betazoom=cn_toepoch(t1,t2,betaEH);
betat=irf_resamp(betazoom,(betazoom(1,1):1/450:betazoom(end,1))');
betat=cn_toepoch(t1,t2,betat);
betatmean=mean(betat,1);
beta=betatmean(2);

c_eval('B?zoom=cn_toepoch(t1,t2,gsmB?);',3:4)
B=mean([B3zoom(:,[1 5]);B4zoom(:,[1 5])]);
B=B(2);

c_eval('Te?zoom=cn_toepoch(t1,t2,eVTe?);',3:4)
Teav=mean(Te4zoom(:,2));

c_eval('Ti?zoom=cn_toepoch(t1,t2,codTi?eV);',4)
Tiav=mean([Ti4zoom(:,2)]);

Vith=irf_plasma_calc(B,Niav,0,Teav,Tiav,'Vtp')/1000; % km/s
Veth=irf_plasma_calc(B,Niav,0,Teav,Tiav,'Vte')/1000; % km/s

r_i=irf_plasma_calc(B,Niav,0,Teav,Tiav,'Rop')/1000; % km
r_e=irf_plasma_calc(B,Niav,0,Teav,Tiav,'Roe')/1000; % km
%%
if 0
d=cn_event_conf(t1,t2,gseB3,gseB4,codVi3,codVi4,gseExB3,gseExB4,gsePos3,gsePos4,eVTe3,codTi4eV,peaNeav,Niav,scpNe3,peaNe3,n_hat,t1str,t2str,'poster');

% Transform into usable quantities when using nhat like basis vector
Bhat=d{1};
M=d{2};
Vi=d{3};
vihat=d{4};
B=d{5};
vdhat=d{6}; % possible propagation direction of wave in field aligned coord
vd=(M\vdhat')'; % possible propagation direction of wave in GSE coord
bh=Bhat;
flh=d{7};
Vith=d{8};
Tiav=d{9};
Teav=d{10};
normdist=d{11};
kdist=d{12};
r_i=d{13};
zdist=d{14};
r_e=d{15};
Ln1=d{16};
end
%%

if 1
    irf_units;
    flh=Units.e*B*1e-9/sqrt(Units.me*Units.mp)/2/pi; % Hz
    window=46;
    flow=0.25;%flh*0.5;
    fs=450;
    [b_av b_hat b_mag]=cn_hat(cn_toepoch(t1,t2,gseB3),cn_toepoch(t1,t2,gseB4));
    z=b_hat;
    y=cn_cross(cn_cross(z,n_hat),z); % close to BL-normal
    x=cn_cross(y,z); % close to drift direction within boundary layer
    gseM=[x;y;z];
    %M=Mamp;
    %gseB3AC=staffB3AC;
    %gseB4AC=staffB4AC;
    c_eval('gseB?AC=irf_filt(cn_toepoch(t1,t2,gseB?),flow,0,fs,5);',3:4)
    %c_eval('gsmB?DC=irf_filt(cn_toepoch(t1,t2,gsmB?),0,flow,fs,5);',3:4)
    c_eval('bB?=cn_m_trans(gseB?AC,gseM,1);',3:4)
    c_eval('gseE?AC=irf_filt(cn_toepoch(t1,t2,gseE?),flow,0,fs,5);',3:4)
    c_eval('bE?=cn_m_trans(cn_toepoch(t1,t2,gseE?),gseM,1);',3:4)
    
    c_eval('bE?AC=cn_m_trans(cn_toepoch(t1,t2,gseE?AC),gseM,1);',3:4)
    
    deltaPos=mean(cn_m_trans(irf_add(1,cn_toepoch(t1,t2,gsePos4),-1,cn_toepoch(t1,t2,gsePos3)),gseM,1));
    deltaPos=deltaPos(2:4);
    
    % Gradient length scale
    Ln2=cn_ln2(t1,t2,t1str,t2str,gsmB3,gsmB4,Teav,Tiav,Neav,Niav,deltaPos(2),deltaPos(1),deltaPos(3),gseM,r_i);
    
    %Ln3=cn_ln2(t1,t2,t1str,t2str,gsmB3DC,gsmB4DC,Teav,Tiav,peaNeav,Niav,normdist,kdist,M,r_i);
    % correlate e-field in propagation direction
    [t0k dtk tplusk tminusk corrk]=cn_mycorr(bE4AC(:,2),bE3AC(:,2),window,1/(bE4AC(2,1)-bE4AC(1,1)),1.03);
    [t02 corrx corr_max]=cn_xcorr(bE4AC,irf_resamp(bE3AC,bE4AC),window,'x');
    [t0kDC dtkDC tpluskDC tminuskDC corrkDC]=cn_mycorr(bE4(:,2),bE3(:,2),window,1/(bE4(2,1)-bE4(1,1)),1.03);
    [t02DC corrxDC corr_maxDC]=cn_xcorr(detrend(bE4,'constant'),detrend(irf_resamp(bE3,bE4),'constant'),window,'x');
    [t02DCy corrxDCy corr_maxDCy]=cn_xcorr(detrend(bE4,'constant'),detrend(irf_resamp(bE3,bE4),'constant'),window,'y');
    t0choice=t02DC;
    
    vchoice=-deltaPos(1)/t0choice;
    
    c_eval('phi?=cn_potential(cn_toepoch(t1,t2,gseE?AC),''GSE'',x*vchoice,''GSE'');',3:4)
    %vpar=zdist/t0choice;
    %c_eval('phi?par=cn_potential(gseE?AC)',3:4)
end
%vchoice=1500;
% Compare dB and potential
Cstr='4';
comp_dB_phi_article3

%Cstr='3';
%comp_dB_phi_article3

propdir='rightleft';
propdir='topbot';
%%
cn_ebplot2(t1,t2,bE3AC,bE4AC,bB3,bB4,...
    gsePos3,gsePos4,vchoice,gseM,propdir,r_e,n_hat,12);