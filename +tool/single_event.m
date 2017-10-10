% Assume data is loaded.
% Also assume the spacecraft is given and that the timeinterval is defined
% in som way beforehand. 

c_eval('n_loc = irf_resamp(peaNe?hf,tint(1),''nearest'');',sc); n_loc = n_loc(2);
di_loc = 9e3*sqrt(n_loc)*sqrt(1/1836); % km
c_eval('re_loc = irf_resamp(re?,tint(1),''nearest'');',sc); re_loc = re_loc(2);
try c_eval('ri_loc = irf_resamp(ri?,tint(1),''nearest'');',sc); ri_loc = ri_loc(2); end
c_eval('flh_loc = irf_resamp(flh?,tint(1),''nearest'');',sc); flh_loc = flh_loc(2);
angles=1:3:360;
f_highpass = flim*flh_loc;
c_eval(['E = irf_tlim(' csys 'E?,tint);'],sc)
c_eval(['B = irf_tlim(' csys 'B?,tint);'],sc);

% Propagation direction
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

avEn=mean(En(:,1+i_dir));
avEk=mean(Ek(:,1+i_dir));
vEn = avEn*y(i_dir,:);
vEk = avEk*x(i_dir,:);
vE = vEn + vEk;
vExB = cross(vE,z(i_dir,:))/B0*1e3;
avExBk=avEn/B0*1e3; % km/s 
avExBn=-avEk/B0*1e3; % km/s 
avExB = sqrt(avExBk^2+avExBn^2);
avExBkdir = x(i_dir,:);
avExBndir = y(i_dir,:);
avExBdir = irf_norm((avExBkdir*avExBk+avExBkdir*avExBn)/irf_abs(avExBkdir*avExBk+avExBndir*avExBn,1));
xx=x(i_dir,:);
yy=y(i_dir,:);
%% Velocity
% Find error
%corr_dir
%withinerror = find(corr_dir>max(corr_dir)*0.95);
Cerror = numel(find(corr_dir>max(corr_dir)*0.95))*3/2;

% Get maximum/minimum variance direction for the eletric field
[~,mva_l,mva_v]=irf_minvar(irf_tlim(E,tint));


% Approximate v range
mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3
v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
v = v_approx*linspace(0.5,1.5,30);
[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
i_v=find(corr_v(1,:)==min(corr_v(1,:)));
velocity=v(i_v);
