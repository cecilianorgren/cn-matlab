function [x,y,z,v_x,v_y,v_z] = f_kepl2svec(GM,a,e,w,W,i,M_0,t_0,t)
%%
% F_KEPL2SVEC computes the state vector (position and velocity) of a
% satellite in a given epoch, from the six keplerian elements. The state 
% vector is refered to the Earth Centered Inertial (ECI) Reference Frame.
%
% HOW: [x,y,z,v_x,v_y,v_z] = f_kepl2svec(GM,a,e,omega,W,i,M_0,t_0,t)
%
% Input:  GM            [1 x 1] product of central body mass and 
%                               gravitational constant in m^3/s^2.
%
%         a             [1 x 1] semi-major axis in meters.
%
%         e             [1 x 1] eccentricity.
%
%         w             [1 x 1] argument of perigee in degrees.
%
%         W             [1 x 1] longitude of ascending node in degrees.
%
%         i             [1 x 1] inclination in degrees.
%
%         M_0           [1 x 1] mean anomaly in degrees at epoch t_0.
%
%         t_0           [1 x 1] epoch of reference in Julian Date.
%
%         t             [n x 1] epoch of computation in Julian Date.
%
% Output: x             [n x 1] x-coordinate in meters.
%
%         y             [n x 1] y-coordinate in meters.
%
%         z             [n x 1] z-coordinate in meters.
%
%         v_x           [n x 1] x-velocity in m/s.
%
%         v_y           [n x 1] y-velocity in m/s.
%
%         v_z           [n x 1] z-velocity in m/s.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 15/12/2015
%
% uses m-files: f_atan2

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 9
    error('Wrong number of input arguments')
end

if max(t_0 > t) == 1
    error('Time epochs for orbit prediction should be greater or equal than the epoch of reference')
end

%Convert degrees to rad
w              = w.*pi/180;
W              = W.*pi/180;
i              = i.*pi/180;
M_0            = M_0.*pi/180;

%Determine the time difference in seconds
dt             = 86400.*(t-t_0);

%Calculate mean anomaly
M              = M_0 + dt.*sqrt(GM/(a^3));

%Normalize mean anomaly
while max(M >= 2*pi); M(M >= 2*pi) = M(M >= 2*pi)-2*pi; end
while max(M < 0)    ; M(M < 0)     = M(M < 0)    +2*pi; end

%Calculate eccentric anomaly using Newton-Raphson
threshold      = 1e-14;

%Initialize eccentric anomaly and difference
E_0            = M;
E_diff         = 1;
ii             = 0;

while max(E_diff) > threshold
    
    E          = E_0 -(E_0-e.*sin(E_0)-M)./(1-e.*cos(E_0));
    E_diff     = abs(E-E_0);
    E_0        = E;
    ii         = ii+1;
    
    if ii == 10000
        
        E_diff = threshold;
        E      = NaN(size(E));
        
    end
    
end

%Calculate true anomaly in rad
v              = 2*f_atan2(sqrt(1 + e).*sin(E/2),sqrt(1 - e).*cos(E/2));

%Calculate eccentric anomaly in meters
r              = a.*(1 - e.*cos(E));

%Calculate position and velocity in Orbital Reference Frame (ORF)
x_orf          = r.*cos(v);
y_orf          = r.*sin(v);
z_orf          = zeros(size(x_orf));

v_x_orf        = (sqrt(GM*a)./r).*(-sin(E));
v_y_orf        = (sqrt(GM*a)./r).*(sqrt(1 - e^2).*cos(E));
v_z_orf        = zeros(size(v_x_orf));

P_orf          = [x_orf,y_orf,z_orf]';
V_orf          = [v_x_orf,v_y_orf,v_z_orf]';

%Set the rotation matrix from ORF to ECI Reference Frame
R(1,1)         = cos(w)*cos(W) - sin(w)*cos(i)*sin(W);
R(1,2)         =-sin(w)*cos(W) - cos(w)*cos(i)*sin(W);
R(1,3)         = sin(i)*sin(W);
R(2,1)         = cos(w)*sin(W) + sin(w)*cos(i)*cos(W);
R(2,2)         =-sin(w)*sin(W) + cos(w)*cos(i)*cos(W);
R(2,3)         =-sin(i)*cos(W);
R(3,1)         = sin(w)*sin(i);
R(3,2)         = cos(w)*sin(i);
R(3,3)         = cos(i);

%Perform the rotation from ORF to ECI Reference Frame
P_eci          = R*P_orf;
V_eci          = R*V_orf;

%Return the results
x              = P_eci(1,:)';
y              = P_eci(2,:)';
z              = P_eci(3,:)';

v_x            = V_eci(1,:)';
v_y            = V_eci(2,:)';
v_z            = V_eci(3,:)';

end
