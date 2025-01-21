function [output] = f_ut12lmst(longitude,julian_day,type)
%%
% F_UT12LMST converts Julian Days to Local Mean Sidereal Time (LMST) for a
% given latitude using the Greenwich Mean Sidereal Time (GMST). More
% information about the algorithm used can be found at: Expressions to
% implement the IAU 2000 definition of UT1, N. Capitaine, P. T. Wallace and
% D. D. McCarthy, A&A, 406, 3, (2003), 1135-1149, DOI:
% http://dx.doi.org/10.1051/0004-6361:20030817.
%
% HOW: [lmst] = f_ut12lmst(longitude,day,month,year,hour,minute,second)
%
% Input:  longitude     [1 x 1] longitude of a place in degrees.
%
%         julian_day    [n x 1] Julian day in decimal days.
%
%         type                  string containing the type of Julian day.
%                               Can be 'JD', 'MJD' or 'TJD'.
%
% Output: output        [n x 1] LMST in degrees.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 16/12/2015
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 3
    error('Wrong number of input arguments')
end

%Calculate Earth's angular velocity
period_Earth                        = (23 + 56/60 + 4.0910/3600)*3600; %sec
omega_Earth                         = 2*pi/period_Earth; %rad/sec

%Convert Julian Day to Gregorian Date
[day,month,year,hour,minute,second] = f_jd2gd(julian_day,type);

%Calculate days elapsed since Julian Day 2451545.0 (1st Jan 2000, 12h UT1)
h_0                                 = zeros(size(day));
d_u                                 = f_gd2jd(day,month,year,h_0,h_0,h_0,'JD')-2451545.0;

%Convert Julian days to Julian centuries of UT1
T_u                                 = d_u/36525;

%Calculate GMST at 0h UT1
gmst_0                              = 67310.54841 + (876600*3600 + 8640184.812866)*T_u + 0.093104*T_u.^2 - 6.2e-6*T_u.^3; %sec

%Correct GMST range
while max(gmst_0) >= 86400; gmst_0(gmst_0 >= 86400) = gmst_0(gmst_0 >= 86400) - 86400; end
while min(gmst_0) < 0     ; gmst_0(gmst_0 < 0)      = gmst_0(gmst_0 < 0) - 86400     ; end

%Convert GMST to degrees
gmst_0                              = gmst_0*(360/24)/3600; %degrees

%Calculate GMST at given time
dt                                  = (f_gd2jd(day,month,year,hour,minute,second,'JD')-f_gd2jd(day,month,year,h_0,h_0,h_0,'JD'))*86400; %sec
gmst                                = gmst_0 + omega_Earth*dt*(180/pi);  %degrees

%Calculate LMST at given time and given longitude
output                              = gmst + longitude; %degrees

%Correct LMST range
while max(output) >= 360; output(output >= 360) = output(output >= 360) - 360; end
while min(output) < 0   ; output(output < 0)    = gmst_0(output < 0) - 360   ; end

end
