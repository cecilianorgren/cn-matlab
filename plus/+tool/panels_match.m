% Load all data
% Needed is E, B, the wanted angles, filtering frequency 

if 1
    cd /Users/Cecilia/Data/BM/20070902/
    load mBS_20070902_1545-1550
    tint=[toepoch([2007 09 02 15 47 29.3]) toepoch([2007 09 02 15 47 30.8])];
else % Load gseE and gseB
    cd /Users/Cecilia/Data/BM/20070831/
    load mBS
    tint = toepoch([2007 08 31 10 19 05.5;2007 08 31 10 19 05.9])';    
end
c_eval('gseB?staff=c_coord_trans(''DSC'',''gse'',dBS?,''cl_id'',?);',3:4);
    c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
    %c_eval('diB?fgm=c_coord_trans(''gse'',''isr2'',gseB?fgm,''cl_id'',?);',3:4);
    c_eval('gseB?fgm=irf_resamp(gseB?fgm,gseB?staff);',3:4);
    c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?staff);',3:4);
    c_eval('gseE?=c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''mat'');',3:4);

%% Do the matching
%t1 = toepoch([2007 08 31 10 19 05.5]);
%t2 = toepoch([2007 08 31 10 19 05.9]);


sc = 3;
c_eval('E = irf_tlim(gseE?,tint);',sc);
c_eval('B = irf_tlim(gseB?,tint);',sc);

angles=1:3:360;
f_highpass=2;
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

%% Velocity and density
n=linspace(0.005,0.6,320); % cc
v=linspace(1,3000,300); % km/s
[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
n_pick=0.4; % cc
n_diff=abs(n-n_pick);
i_n=find(n_diff==min(n_diff));
i_v=find(abs(corr_v(i_n,:))==min(abs(corr_v(i_n,:))));
velocity=v(i_v);

%% Visualize the matching
% Direction
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En); 
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);

%% Velocity
% if more than one densitiy, choose one by specifying index
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v_2.gif','DelayTime',0.01,'LoopCount',inf);

%% Density vs velocity
figure; h=axes;
axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);

%% Direction with panels
nPanels = 12;
axis_handle = irf_match_phibe_vis('direction_panels',x,y,z,corr_dir,intEdt,Bz,nPanels); 

%% Plot potential with length scale
h=axes;
irf_plot(h,phi_E(:,[1 1+i_v]))
ylabel(h,'\phi [V]')
grid(h,'off')
    v1=velocity;
    
    %
    % Left panels
    %
    ax1=h(1);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    tick0=1;
    % put each 2rd tick
    tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    tticks = xticks0(ticks);
    lticks = round(tticks*v1-1*xticks0(1)*v1);
    %lticks = round(tticks*v1-xticks0(tick0)*v1);
    
    %nticksl=xlim(1)/5/r_e;
    %nticksr=xlim(2)/5/r_e;
    %xticks=0:xlim(2)/nticksr:xlim(2);
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end    
    %xticklabels{1} = ' '; 
    %xticklabels{2} = ' ';
    xticklabels{end} = ' '; 
    xticklabels{end} = 'km'; 
    yticks=[];
    yticklabels=[];
    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
             'XAxisLocation','top','YAxisLocation','right',...
             'Color','none','Ytick',yticks,...
             'YTickLabel',yticklabels,...
             'xlim',xlim,'xtick',tticks,...
             'XTickLabel',xticklabels);

         hold(h,'off')