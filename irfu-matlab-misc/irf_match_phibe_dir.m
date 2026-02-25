function [x y z correlation intEdt Bz B0 dEk dEn Ek En] = irf_match_phibe_dir(B,E,angles,f)
% IRF_MATCH_PHIBE_DIR Get propagation direction by matching dBpar and "phi".
%   Tries different propagation directions and finds the direction 
%   perpendicular to the magnetic field that gives the best correlation
%   between the electrostatic potential and the parallel wave magnetic field
%   according to int(E)dt = waveB*B0/n*e*mu0.
%
%   [x,y,z,correlation,intEdt,Bz,dEk,dEn,Ek,En] = IRF_MATCH_PHIBE_DIR(B,E,angles,f)
%
%   Input 
%       B - magnetic field (to be filtered if f is given)
%           can be given in cell where B = {BDC,BAC} or {fluxgate B,searchcoil B} or {[B0x,B0y,B0z],BAC}
%       E - electric field (to be filtered if f is given)
%       angles - the angles in degrees to try (1-180 default)
%       f - filter frequency
% 
%   Output 
%       x - normal direction (size: n_tries x 3)
%       y - propagation direction
%       z - magnetic field direction.
%       correlation - correlation vector
%       intEdt - "potential"
%       Bz - wave magnetic field in parallel direction
%       B0 - mean magnetic field
%       dEk - wave electric field in propagation direction
%       dEn - wave electric field in propagation normal direction
%       Ek - electric field in propagation direction
%       En - electric field in propagation normal direction
% 
%   Examples:
%       % Direction
%       angles=1:3:360;
%       f_highpass=7;
%       [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=IRF_MATCH_PHIBE_DIR(B,E,angles,f_highpass);
%       i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
%       direction=x(i_dir,:);
%
%       % Velocity and density
%       n=linspace(0.01,0.1,100);
%       v=linspace(200,2000,100);
%       [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
%       i_v=find(corr_v(:,1)==min(corr_v(:,1)));
%       velocity=v(i_v);
%   
%       % Figures
%       gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,En,Ek);      
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);
%
%       i_n=50; % if more than one densitiy, choose one by specifying index
%       gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
%       imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);
%
%       figure; h=axes;
%       axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
%
%   See also IRF_MATCH_PHIBE_V, IRF_MATCH_PHIBE_VIS

% Check to see if B is a structure, if yes, then it is {B0,waveB or scB}
if iscell(B)    
    BDC = B{1}; % could be fluxgate data
    BAC = B{2}; % could be searchcoil data
else % data can be combined already
    BDC = B;
    BAC = B;
end
if isa(BDC(1),'TSeries'); tsBDC = BDC; BDC = [tsBDC.time.epochUnix double(BDC.data)]; end
if isa(BAC(1),'TSeries'); tsBAC = BAC; BAC = [tsBAC.time.epochUnix double(BAC.data)]; end
if isa(E,'TSeries'); tsE = E; E = [tsE.time.epochUnix double(E.data)]; end
    
% Filter if f is given, otherwise assume it is filtered
if exist('f','var') 
    if numel(f)==1 % highpass filter at given frequency
      fs = 1/(BAC(2,1)-BAC(1,1));
      BAC=irf_filt(BAC,f,0,fs,5);
      fs = 1/(E(2,1)-E(1,1));
      EAC=irf_filt(E,f,0,fs,5);
    else % bandpass filter in given frequency interval
      fs = 1/(BAC(2,1)-BAC(1,1));
      BAC=irf_filt(BAC,f(1),f(2),fs,5);
      fs = 1/(E(2,1)-E(1,1));
      EAC=irf_filt(E,f(1),f(2),fs,5);
    end
else
    EAC=E;
end

% Resample B to E if they have different size
if size(BAC,1) ~= size(EAC,1)
    BAC=irf_resamp(BAC,EAC);
end

% Get background magnetic field, for field aligned direction and irf_match_phibe_v
if size(BDC,2)==3; B0 = norm(BDC); B0hat = irf_norm(mean(BDC(:,1:3),1)); % [Bx,By,Bz]
elseif size(BDC,2)==4; B0=mean(irf_abs(BDC,1)); B0hat = irf_norm(mean(BDC(:,2:4),1)); % [t Bx By Bz]
elseif size(BDC,2)==5; B0=mean(BDC(:,end));  B0hat = irf_norm(mean(BDC(:,2:4),1)); % [t Bx By Bz |B|] - B contains absolute value
end

% If no angles are specified, set 1,4,7,...,158 as default
if ~exist('angles','var') 
    angles=1:3:360;
end
na=length(angles); % number of angles

% Set up coordinate systems
if 1
    z=repmat(B0hat,na,1); % B/z direction, tries*3
    y=irf_norm(irf_cross(irf_cross(z,[1 0 0]),z)); % perp1
    x=irf_norm(irf_cross(y,z)); % perp2
    
    theta=(0:2*pi/na:2*pi-pi/na)'; % angles
    xn=irf_norm(x.*repmat(cos(theta),1,3)+y.*repmat(sin(theta),1,3));
    y=irf_cross(z,xn);
    x=xn;
end
    
% Field aligned B
Bz=irf_dot(BAC,z(1,:));

% Allocate correlations
correlation=zeros(na,1);

% Allocate vectors, 4 first used for illustration
dEk=[EAC(:,1) zeros(size(E,1),na)]; % field to integrate
dEn=[EAC(:,1) zeros(size(E,1),na)]; % normal field
Ek=[E(:,1) zeros(size(E,1),na)]; % unfiltered
En=[E(:,1) zeros(size(E,1),na)]; % unfiltered
intEdt=[EAC(:,1) zeros(size(E,1),na)]; % potential

% Integrate E in all x-directions
for k=1:na
    % Used for visualization
    dEk(:,k+1)=irf_dot(EAC(:,(2:4)),x(k,:));
    dEn(:,k+1)=irf_dot(EAC(:,(2:4)),y(k,:));    
    Ek(:,k+1)=irf_dot(E(:,(2:4)),x(k,:));
    En(:,k+1)=irf_dot(E(:,(2:4)),y(k,:));
    
    % Get Phi_E = int(Ek), there's no minus since the field is integrated
    % in the opposite direction of the wave propagation direction.
    prel=integ([dEk(:,1) dEk(:,k+1)]);    
    intEdt(:,k+1)=prel(:,2)-nanmean(prel(:,2));
       
    % Get correlation
    correlation(k,1)=xcorr(intEdt(:,k+1),Bz(:,2),0,'coeff');    
end

if 0
vis.type='direction';
vis.x=x;
vis.y=y;
vis.z=z;
vis.correlation=correlation;
vis.intEdt=intEdt;
vis.Ek=Ek;
vis.En=En;
vis.ufEk=ufEk;
vis.ufEn=ufEn;
end
end
function xint=integ(x,tref,time_step)
%IRF_INTEGRATE  Integrate time series
%
% xint=irf_integrate(x,tref,time_step)
%   integrate time series. time steps that are larger
%   than 3times the time step are assumed to be data gaps.
%
%   x - time series  to integrate
%   tref - optional, integration starting time (optional) 
%        isdat epoch or [yyyy mm dd hh mm ss.ss]
%   time_step - optional, all time_steps larger than 3*time_step are
%   assumed data gaps, default is that time_step is the smallest value of
%   all time_steps of the time series

dt=[0 ; diff(x(:,1))];
if nargin < 3, % estimate time step
    time_steps=diff(x(:,1));
    [~,ind_min]=min(time_steps);
    time_steps(ind_min)=[]; % remove the smallest time step in case some problems
    time_step=min(time_steps);
end
dt(dt>3*time_step)=0;
xint=x;
for j=2:size(xint,2),
  j_ok=find(~isnan(xint(:,j)));
  %xint(j_ok,j)=cumsum(x(j_ok,j).*dt(j_ok),1);
  xint(j_ok,j)=cumtrapz(x(j_ok,j).*dt(j_ok),1);  
end

if nargin>=2, % other origo for integration 
    if size(tref)==6,
        tt=toepoch(tref);
    elseif size(tref)==1;
        tt=tref;
    else
        errS = 'do not know how to treat TREF input';
        irf.log('critical',errS), error(errS)
    end
    if tt < x(1,1),tt=x(1,1);end % if tref before first data point,set it to time of first data point
    if tt > x(end,1),tt=x(end,1);end % if tref after laast data point,set it to time of last data point
    xint_ref=irf_resamp(xint,tt,'nearest');
    xint=irf_add(1,xint,-1,xint_ref);
end
end