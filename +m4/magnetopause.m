filepath = '/Users/Cecilia/M4/thor/orbit_coverage/';
fid = fopen([filepath 'bowshock_magnetopause_pos_solarwind.txt']);
format ='%*s%*s%f%f%*f';
A = textscan(fid,format);
R_BS = sort(A{1}); % RE
R_MP = sort(A{2}); % RE
%V_SW = A{3}; % km/s
N = numel(R_BS);
indR = round([1 0.25*N 0.5*N 0.75*N N]);
linestyle = {'--','--','-','--','--'};
%%
ii=3;
R0 = R_BS(indR(ii));
% Magnetopause
Bz = 0; Dp = 2;
theta=0:0.1:pi;
alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
r=R0*(2./(1+cos(theta))).^alpha;
xMP=r.*cos(theta);
yMP=r.*sin(theta);
yMP(abs(xMP)>100)=[];
xMP(abs(xMP)>100)=[];

% Bowshock
xBS=R0:-0.5:-100;
yBS=sqrt(0.04*(xBS-R0).^2-45.3*(xBS-R0)); % original F/G model adds rstandoff^2=645
