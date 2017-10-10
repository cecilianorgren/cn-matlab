function [x0,y0,x1,y1,m0,m2] = newposition(vx0,vy0,initialization,rlim,B0)
% Randomizes a wighted radial position and shifts it with the gyroradius so
% that it is in fact the gyrocenter that was decided. 
%
%   [x1,y1,x2,y2,m1,m2] = newposition(N,vx0,vy0);       
%       vx, vy - vectors with their thermal velocity given in km/s
%       x1, y1 - directly generated positions depending on model
%       x2, y2 - position a gyroradius shifted so that x0, y0 besomces the gyrocenter
%       m1     - the particle weight, depending on r0 = sqrt(x0^2+y0^2)
%       m2     - the particle weight, depending on r1 = sqrt(x1^2+y1^2)

N = numel(vx0);

switch initialization
    case 'normal'
        r0 = rlim*randn(N,1);
        angle = 180*rand(N,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);        
        m0 = 2*pi*rlim*exp(r0.^2/rlim.^2/2);
        disp('Using normal distribution for starting positions.')
    case 'r'
        r0 = rlim*rand(N,1);
        angle = 360*rand(N,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);        
        m0 = r0/rlim; 
        disp('Using r distribution for starting positions.')       
    case 'rm'
        r0 = rlim*rand(N,1);
        angle = 360*rand(N,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);        
        m0 = r0/rlim; 
        disp('Using r distribution for starting positions.')       
    case 'uniform'
        % pdf = r;
        % cdf = r^2/2
        cdfy = rand(N,1);
        r0 = 2*sqrt(cdfy);
        r0 = r0/max(r0)*rlim;
        %r0 = rlim*rand(N,1);
        angle = 360*rand(N,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);  
        r0 = sqrt(x0.^2+y0.^2);
        m0 = 0.1*ones(N,1); 
        disp('Using uniform distribution for starting positions.')
end

vperp = sqrt(vx0.^2+vy0.^2); % km/s
me = 9.10939999999e-31;
e = 1.6022e-19;
rg0 = me*vperp/e/(B0*1e-9); % km, the starting position should be shifted 
                            % with this distance
vxg0 = find(vx0>=0); % right half xy-plane
vxs0 = find(vx0<0); % left half xy-plane
vth = zeros(size(vx0));
vth(vxg0) = atand(vy0(vxg0)./vx0(vxg0));
vth(vxs0) = atand(vy0(vxs0)./vx0(vxs0)) + 180; 
vth = vth-90;

% Shift thestarting position from the gyrocenter
x1 = x0 + rg0.*cosd(vth);
y1 = y0 + rg0.*sind(vth);

r1 = sqrt(x0.^2+y0.^2);
switch initialization
    case 'normal'              
        m2 = 2*pi*rlim*exp(r1.^2/rlim.^2/2);
        %disp('Using normal distribution for starting positions.')
    case 'r'
        m2 = r1/rlim; 
    case 'rm'
        m2 = r0/rlim;         
        %disp('Using r distribution for starting positions.')       
    case 'uniform'              
        m2 = 0.1*ones(N,1);         
end
