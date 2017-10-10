

% load electric field
cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields

if ~exist('ind','var'); ind = 5; end
ind = 5; % Has nice mono.
tint = tint{highQuality(ind)};

epar4 = irf_tlim(EparAC4,tint);
epar3 = irf_tlim(EparAC3,tint);
eper4 = irf_tlim(Epar4,tint);
eper3 = irf_tlim(Epar3,tint);

dx=cn.mean(irf_tlim(dp1,tint),1); % SP
dy=cn.mean(irf_tlim(dp2,tint),1); % BxSP
dz=cn.mean(irf_tlim(dp3,tint),1); % B


% model from Chen2004: Bernstein-Greene-Kruskal solitary waves in...
phi = @(phi0,r,z,lz,lr) phi0*exp(-z.^2/2/lz.^2-r.^2/2/lr.^2);
lz = 10; % km
lr = 10; % km
phi0 = 300; % V

% try different sc postitions
% R4 follows from R3
x3 = 10; x4 = x3-dx; % in spin plane
y3 = 10; y4 = y3-dy; % out of spin plane

%% make plot
if 1 % sc position
    % make plot
    x=linspace(-lr*0.60,2*lr,200);
    y=linspace(-2*lr,2*lr,200);
    [X,Y]=meshgrid(x,y);

    R = sqrt(X.^2+Y.^2);
    levels = phi0*[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 ];
    [C,h] = contour(X,Y,phi(phi0,R,0,lz,lr),levels);
    clabel(C,h,levels(1:2:end),'labelspacing',1000);
    xlabel('in sp (\perp B)')
    ylabel('out of sp (\perp B)')
    hold on;
    plot(x3,y3,'og',x4,y4,'vb','markersize',20,'linewidth',2)
end
if 1 % fit electric field

zlim = 30; % km
nz = 100;
z = linspace(-zlim,zlim,nz);
theta = linspace(0,359,360);
rlim = 30; nr = 50;
r = linspace(0,rlim,nr);
x = r'*cosd(theta);
y = r'*sind(theta);

[X Y Z] = meshgrid(x,y,z);
% define trajectory, radius at which eh pass
rt = 5; % km