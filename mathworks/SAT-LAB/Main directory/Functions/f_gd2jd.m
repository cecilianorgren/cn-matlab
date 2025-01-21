function [output] = f_gd2jd(day,month,year,hour,minute,second,type)
%%
% F_GD2JD converts Gregorian Dates to Julian Days of many formats. The UT 
% hour is also used to calculate the decimal part of the day, which is
% refered to Mean Solar Time. More information about the algorithm used can
% be found at http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html .
%
% HOW: [jd] = f_gd2jd(day,month,year,hours,minutes,seconds,type)
%
% Input:  day           [n x 1] day of Gregorian Calendar.
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
%         type                  string containing the type of Julian day.
%                               Can be 'JD', 'MJD' or 'TJD'.
%
% Output: output        [n x 1] Julian day in decimal days.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 16/12/2015
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 7
    error('Wrong number of input arguments')
end

if strcmp(type,'JD') == 0 && strcmp(type,'MJD') == 0 && strcmp(type,'TJD') == 0
    error('type parameter should be "JD", "MJD" or "TJD"')
end

%If the month is January or February, subtract 1 from the year and add 12 
%to the month.
year(month<3)  = year(month<3) - 1;
month(month<3) = month(month<3) + 12;

%Force the elements to have a collumn vector form
year           = year(:);
month          = month(:);

%Calculate intermediate variables
A              = floor(year/100);
B              = floor(A/4);
C              = 2 - A + B;
D              = floor(365.25*(year + 4716));
E              = floor(30.6001*(month + 1));
jdn            = C + D + E + day - 1524.5;

%Calculate the decimal part of Julian Day
output         = jdn + hour/24 + minute/1440 + second/86400;

if strcmp(type,'MJD') == 1
    
    %Convert to Modified Julian Day
    output     = output - 2400000.5;
    
elseif strcmp(type,'TJD') == 1
    
    %Convert to Truncated Julian Day
    output     = floor(output - 2440000.5);
    
end

end
