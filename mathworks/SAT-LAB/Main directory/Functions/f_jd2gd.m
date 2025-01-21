function [day,month,year,hour,minute,second] = f_jd2gd(julian_day,type)
%%
% F_JD2GD converts Julian Days of many formats to Gregorian Dates. More
% information about the algorithm used can be found at:
% http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html .
%
% HOW: [day,month,year,hour,minute,second] = f_jd2gd(julian_day,type)
%
% Input:  julian_day    [n x 1] Julian day in decimal days.
%
%         type                  string containing the type of Julian day.
%                               Can be 'JD', 'MJD' or 'TJD'.
%
% Output: day           [n x 1] day of Gregorian Calendar.
%
%         month         [n x 1] month of Gregorian Calendar.
%
%         year          [n x 1] year of Gregorian Calendar.
%
%         hour          [n x 1] hour part of UT.
%
%         minute        [n x 1] minute part of UT.
%
%         second        [n x 1] second part of UT.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 18/12/2015
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 2
    error('Wrong number of input arguments')
end

if strcmp(type,'JD') == 0 && strcmp(type,'MJD') == 0 && strcmp(type,'TJD') == 0
    error('type parameter should be "JD", "MJD" or "TJD"')
end

%Convert input day to Julian Day
if strcmp(type,'MJD') == 1
    
    %Convert Modified Julian Day to Julian Day
    julian_day   = julian_day + 2400000.5;
    
elseif strcmp(type,'TJD') == 1
    
    %Convert Truncated Julian Day to Julian Day
    julian_day   = julian_day + 2440000.5;
    
end

%Calculate intermediate variables
Q                = julian_day + 0.5;
Z                = floor(Q);
W                = floor((Z - 1867216.25)/36524.25);
X                = floor(W/4);
A                = Z + 1 + W - X;
B                = A + 1524;
C                = floor((B - 122.1)/365.25);
D                = floor(365.25*C);
E                = floor((B - D)/30.6001);
F                = floor(30.6001*E);

%Calculate day, month and year
day              = B - D - F + (Q - Z);
month            = E - 13;

%Correct month range
while max(month) > 12; month(month > 12) = month(month > 12) - 12 ; end
while min(month) < 1 ; month(month < 1)  = month(month < 1) + 12  ; end

year(month < 3)  = C(month < 3) - 4715;
year(month >= 3) = C(month >= 3) - 4716;

%Calculate hour, minute and second
hour             = (Q - Z)*24;
minute           = (hour - floor(hour))*60;
second           = (minute - floor(minute))*60;

%Get the integer part of day, hour and minute
day              = floor(day);
hour             = floor(hour);
minute           = floor(minute);

end
