function [A,im,map]=irf_plasma_wave_visualization(varargin)
% IRF_PLASMA_WAVE_VISUALIZATION visualize cold plasma waves
% [A,im,map]=IRF_PLASMA_WAVE_VISUALIZATION(...)
%
% IRF_PLASMA_WAVE_VISUALIZATION(Parameter,ParameterValue,...) - specify input
%
% IRF_PLASMA_WAVE_VISUALIZATION(AX,...) - visualize in specified axis
%
% Reference system: Z is along field ambient field B0
%                   k=[kx 0 kz], wave vector has 0 Y-component
%
% Input is pairs: 'ParameterName', ParameterValue
%  Available parameters
%       'E' - wave polarization  [Ex Ey Ez] or [Ex Ez] or [Ex]
%      'kx' - X component of wave vector (default 2*pi, 1 wave period)
%      'kz' - Z component of wave vector (default 0)
%      'qm' - vector of charge/mass for each species (default [1 -1836])
%       'N' - number of particles in each species
%       'f' - wave frequency in Hz (default 1, 20 frames per s)
%       'T' - length of simulation
%       'X' - X length (default is 1)
%       'Z' - Z length (default is 1)
%'init_type'- 0 - fixed square, 1 - extended square, 2- zoomed square
%  'Nfield' - number of field lines (default 20)
%'Nfieldlinepoints' - number of points to resolve each field line (default 20)
%    'view' - 'XZ','XY','YZ','3D' or vector as defined by matlab get(gca,'view')
%    'demo' - available: 'alfven_compr','alfven_compr45','alfven_shear','whistler','alfven_shear'
%
% Output: 	A - movie matrix
%           im - image matrix
%           map - color map
%
% Examples:
%   [A,im,map]=irf_plasma_wave_visualization('demo','alfven_shear');
%   movie(A,20);
%
% to create animated gif
% imwrite(im,map,'anim.gif','DelayTime',0,'LoopCount',inf)

% For 3D visulatization one needs to introduce E,k,q/m instead of B, k, v.
% Because for example electrostatic wave will also have associated field line
% motion which one cannot obtain from k and B

global Xplot_lim Zplot_lim flag_plot_initiated hlines hparticles
flag_plot_initiated=0;
flag_axis_view_given=0;
%% check output
flag_generate_output_files=0;
if nargout>=1, flag_generate_output_files=1;end
%% default values that can be overwritten by input
B0 = 1; % default background magnetic field is 1
init_type=0;
X=1;Y=1;Z=1;

axis_view=[-45 45]; % 'XZ' view

% overall default parameters
N=1000;
Nfield=20; % number of field lines
Nfieldlinepoints=20; % number of field line points

% default upstream parameters
Bn1=1; % default background magnetic field upstream discontinuity
Bt1=1;
vn1=0; % relative velocity of discontinuity and particles
vt1=1;
p1=1;
rho1=1;
rhovn1=1;

% default jump in parameters
dBn=0;
dBt=0; % difference in relative tangential magnetic field
dvn=0;
dvt=0; % difference in relative tangential velocity
dp=0; % difference in pressure
drho=1; % difference in density
drhovn=0;

% default downstream parameters
Bn2=Bn1+dBn;
Bt2=Bt1+dBt;
vn2=vn1+dvn;
vt2=vt1+dvt;
p1=p1+dp;
tho2=rho1+drho;
rhovn2=rhovn1+drhovn;

%% check input axes
[ax,args,nargs] = axescheck(varargin{:});
if 0 %% check input
if nargs == 0 , % show help
    help irf_plasma_wave_visualization;
    return
end
if (rem(nargs,2) ~= 0), % input not in pair forms, demos if necessary
    error('IRFU_MATLAB:irf_plasma_wave_visualization:InvalidNumberOfInputs','Incorrect number of input arguments')
end
while 1
    if isempty(args), break, end
    l = 2;
    switch lower(args{1})
        case 'z'
            if isnumeric(args{2})
                Z = args{2};
            else irf_log('fcal,','wrongArgType : Z must be numeric')
            end
        case 'x'
            if isnumeric(args{2})
                X = args{2};
            else irf_log('fcal,','wrongArgType : Z must be numeric')
            end
        case 't'
            if isnumeric(args{2})
                T = args{2};
            else irf_log('fcal,','wrongArgType : T must be numeric')
            end
        case 'n' % number of particles in every species
            if isnumeric(args{2})
                N = args{2};
            else irf_log('fcal,','wrongArgType : N must be numeric')
            end             
        case 'nfield' % number of field lines
            if isnumeric(args{2})
                Nfield = args{2};
            else irf_log('fcal,','wrongArgType : Nfield must be numeric')
            end
        case 'nfieldlinepoints' % number of field lines
            if isnumeric(args{2})
                Nfieldlinepoints = args{2};
            else irf_log('fcal,','wrongArgType : Nfieldlinepoints must be numeric')
            end
        case 'init_type'
            if isnumeric(args{2})
                init_type = args{2};
            else irf_log('fcal,','wrongArgType : init_type must be numeric')
            end
        case 'view'
            if ischar(args{2})
                switch lower(args{2}), % expecting 'XZ','XY' or 'YZ'
                    case 'xy'
                        axis_view=[0 90];
                        flag_axis_view_given=1;
                    case 'yz'
                        axis_view=[90 0];
                        flag_axis_view_given=1;
                    case 'xz'
                        axis_view=[0 0];
                        flag_axis_view_given=1;
                    case '3d'
                        axis_view=[-45 45];
                        flag_axis_view_given=1;
                    otherwise
                        irf_log('fcal,','wrongArgType : view string can be XY, XZ or YZ')
                end
            elseif isnumeric(args{2}), % specifies AZ EL
                axis_view=args{2};
                flag_axis_view_given=1;
            else
                irf_log('fcal,','wrongArgType : ''view'' has wrong parameter value')
            end
        case 'demo'
            if ischar(args{2})
                switch lower(args{2}), % expecting 'XZ','XY' or 'YZ'
                    case 'contact_discontinuity'
                    case 'tangential_discontinuity'
                    case 'rotational_discontinuity' 
                    case 'intermediate_shock'  
                    case 'fast_shock'
                    case 'slow_shock'
                    otherwise
                        irf_log('fcal,','wrongArgType : ''demo'' parameter value unknown')
                end
            else
                irf_log('fcal,','wrongArgType : ''demo'' has wrong parameter value')
            end
    end
    args = args(l+1:end);
end
end
%% initialize figure
if isempty(ax), % axis not defined
    set(0,'defaultLineLineWidth', 1.5);
    if any(get(0,'children') == 101), % reuse existing window
        figure(101);
        cla;
        h=gca;
        set(h,'view',axis_view);
    else                              % create new window
        fn=figure(101);clf;
        set(fn,'Renderer','painters')
        set(fn,'color','white')
        set(gcf,'PaperUnits','centimeters')
        xSize = 12; ySize = 12;
        xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        pos=get(gcf,'Position');
        set(gcf,'Position',[pos(1) pos(2) xSize*60 ySize*60])
        h=axes('position',[0 0 1 1]); % [x y dx dy]
        cla;
        set(h,'view',axis_view);
        set(h,'visible','off')
    end
else % use speciifed axis
    cla;
    h=gca;
end
%% construct dB,V,dR disturbances and gyrofrequencies

% construct magnetic field disturbance as discontinuity passes
% 

V=zeros(numel(qm),3);     % allocate particle velocity disturbance
for j=1:numel(qm), % step through all particles
    Omega_c(j)=qm(j)*B0; %  gyrofrequencies
    V(j,1) = 1i* qm(j) / (2*pi*f.^2 - Omega_c(j)^2) ...
        * (2*pi*f*E(1) + 1i*Omega_c(j)*E(2)); % x component
    V(j,2) = 1i* qm(j) / (2*pi*f.^2 - Omega_c(j)^2) ...
        * (-1i*Omega_c(j)*E(1) + 2*pi*f*E(2)); % y component
    V(j,3) = 1i* qm(j) / (2*pi*f) * E(3); % z component
end
dR=1i*V/(2*pi*f); % position disturbances

% construct field line disturbance
% w = (E x B) / B^2
% w = - i 2 pi f r
% r = i / (2 pi f) w =  i / (2 pi f) (ExB) / B^2
bvec= [0 0 B0]/B0^2;
dRfieldline=1i*cross(E,bvec)/(2*pi*f); % field line disturbance vector
for j=1:numel(qm), % display information on disturbances
    fprintf('%d. ===== q/m = %5.1f =====\n',j,qm(j));
    disp(['dR = ' num2str(dR(j,:),'%6.2f')]);
    disp(['dV = ' num2str(V(j,:),'%6.2f')]);
end
disp('3. ===== field lines =====');
disp(['dR = ' num2str(dRfieldline,'%6.2f')]);
%% initialize location of particles
r_init=cell(numel(qm),1); % initial positions cell for each species
if init_type==0 || init_type==1, % random particles in
    for j=1:numel(qm)
        r_init{j}=rand(N,3); % random positions
        r_init{j}(:,1)=r_init{j}(:,1)*X;
        r_init{j}(:,3)=r_init{j}(:,3)*Z;
        r_init{j}(:,2)=r_init{j}(:,2)*0;
        %        dVxlim=max(abs(dVi(1)),abs(dVi(1)));
        %        dVylim=max(abs(dVi(2)),abs(dVi(2)));
    end
    if init_type==0,
        Xplot_lim=[0 X];
        Yplot_lim=[-Y/2 Y/2];
        Zplot_lim=[0 Z];
    elseif init_type==1,
        dVxlim=max(abs(V(:,1)));
        dVylim=max(abs(V(:,2)));
        dVzlim=max(abs(V(:,3)));
        Xplot_lim=[-dVxlim/f X+dVxlim/f];
        Yplot_lim=[-Y/2-dVylim/f Y/2+dVylim/f];
        Zplot_lim=[-dVzlim/f Z+dVzlim/f];
    end
elseif init_type==2, % spread depending on disturbance amplitude
    for j=1:numel(qm)
        r_init{j}=rand(N,3); % random positions
        r_init{j}(:,1)=r_init{j}(:,1)*(X+2*abs(V(j,1))/real(f))-abs(V(j,1))/real(f);
        r_init{j}(:,3)=r_init{j}(:,3)*(Z+2*abs(V(j,3))/real(f))-abs(V(j,3))/real(f);
        r_init{j}(:,2)=r_init{j}(:,2)*0;
        %        dVxlim=max(abs(dVi(1)),abs(dVi(1)));
        %        dVylim=max(abs(dVi(2)),abs(dVi(2)));
    end
    Xplot_lim=[0 X];
    Yplot_lim=[-Y/2 Y/2];
    Zplot_lim=[0 Z];
end
%% initialize field lines
field_lines_init=zeros(Nfield,Nfieldlinepoints,3); % x and y coordinates of field lines
xfield=0:X/Nfield:X;
for j=1:Nfield,
    field_lines_init(j,:,1)=xfield(j);
    field_lines_init(j,:,2)=0;
    field_lines_init(j,:,3)=0:Z/(Nfieldlinepoints-1):Z;
end
%% plot initial status
% showinitial statuts
plot_particles(h,r_init,qm,field_lines_init);
set(h,'xlim',Xplot_lim,'ylim',Yplot_lim,'zlim',Zplot_lim);
%xlabel(h,'x');ylabel(h,'y');zlabel(h,'z');
hold(h,'on');
flag_plot_initiated=1;
if ~flag_axis_view_given, % allow to adjust axes interactively
    input('Press enter to continue:','s');
end
%% step through wave
dt=(1/abs(f))/30; % 20 frames per period
time_steps=0:dt:T;
r=cell(length(time_steps),1); % initialize position of particles
for jj=1:numel(qm),
    k_phase_factor=1i*dot(repmat(k,N,1),r_init{jj},2);
    rep_dR=repmat(dR(jj,:),N,1);
    for j=1:length(time_steps), % calculate particle positions
        t=time_steps(j);
        phasefactor=exp(-1i*2*pi*f*t+k_phase_factor); % size Npart x 1
        rep_phasefactor=repmat(phasefactor,1,3); % size Npart x 3
        r{j}{jj}=r_init{jj} + real(rep_dR.*rep_phasefactor);
    end
end
if 1, % calculate field line position
    % v = -i w r > r = i v / w
    % v = E x B / B^2
    % B x v = [B^2 E - (E.B)B] / B^2
    % k x E = -w dB
    % (k.E) k - k^2 E = -w dB x k
    
    field_lines=zeros(length(time_steps),Nfield,Nfieldlinepoints,3); % field line time evolution
    rep_k=repmat(shiftdim(k,-1),[Nfield,Nfieldlinepoints,1]); % replicated k vector the size of field_lines_init
    rep_dRfieldline = repmat((shiftdim(repmat(dRfieldline,Nfieldlinepoints,1),-1)),[Nfield 1 1]); % replicated to size Nfield x Nfieldlinepoints x 3
    k_phase_factor=1i*dot(rep_k,field_lines_init,3); % position phase factor
    for j=1:length(time_steps), % calculate field lines position
        t=time_steps(j);
        phasefactor = exp(-1i*2*pi*f*t+k_phase_factor); % size Nfield x Nfieldlinepoints
        rep_phasefactor = shiftdim(repmat(shiftdim(phasefactor,-1),3,1),1); % create phase in matrix Nfield x Nfieldlinepoints x 3
        dr_field=real(rep_dRfieldline.*rep_phasefactor);
        field_lines(j,:,:,:)=field_lines_init+dr_field;
    end
end
if 1, % plot wave animation
    A=moviein(length(time_steps)); % create matlab movie matrix
    for j=1:length(time_steps), % plot wave animation
        %pause(.01);
        plot_particles(h,r{j},qm,squeeze(field_lines(j,:,:,:)));
        drawnow;
        if flag_generate_output_files,
            f=getframe(h);
            A(:,j)=f;
            if j==1, % initialize animated gif matrix
                [im,map] = rgb2ind(f.cdata,256,'nodither');
                im(1,1,1,20) = 0;
            else
                im(:,:,1,j) = rgb2ind(f.cdata,map,'nodither');
            end
        end
    end
end

    function  plot_particles(h,r,qm,field_lines,markers)
        persistent marks
        if nargin==5, marks=markers; end % marks passed as argument
        if nargin<4, field_lines=[]; end
        if nargin<1, h=gca;end
        if ~flag_plot_initiated,      % first run
            if nargin<5,              % marks not defined
                for ii=1:numel(qm),
                    marks{ii}.markersize=10+10*((log(qm(ii))-log(min(qm)))/(log(max(qm))-log(min(qm))));
                    if sign(qm(ii))>0 % positive charge
                        marks{ii}.color='red';
                    else % negative charge
                        marks{ii}.color='blue';
                    end
                    marks{ii}.marker='.';
                end
            end
            cla(h);
            for ii=1:length(r),       % plot particles
                hparticles(ii)=line(r{ii}(:,1),r{ii}(:,2),r{ii}(:,3),'linestyle','none',...
                    'markersize',marks{ii}.markersize,...
                    'marker',marks{ii}.marker,...
                    'color',marks{ii}.color);
            end
            if ~isempty(field_lines), % plot field lines
                for jjj=1:size(field_lines,1),
                    hlines(jjj)=line(field_lines(jjj,:,1),field_lines(jjj,:,2),field_lines(jjj,:,3),'linestyle','-',...
                        'color','black');
                end
            end
        else                          % repeated run, reuse handles
            for ii=1:length(r),       % replot particles
                delete(hparticles(ii));
                hparticles(ii)=line(r{ii}(:,1),r{ii}(:,2),r{ii}(:,3),'linestyle','none',...
                    'markersize',marks{ii}.markersize,...
                    'marker',marks{ii}.marker,...
                    'color',marks{ii}.color);
            end
            if ~isempty(field_lines), % replot field lines
                for jjj=1:size(field_lines,1),
                    delete(hlines(jjj));
                    hlines(jjj)=line(field_lines(jjj,:,1),field_lines(jjj,:,2),field_lines(jjj,:,3),'linestyle','-',...
                        'color','black');
                end
            end
        end
    end
end
