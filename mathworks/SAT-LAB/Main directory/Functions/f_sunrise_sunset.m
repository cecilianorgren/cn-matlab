function [UT_sunrise, UT_sunset] = f_sunrise_sunset(day_of_year,latitude,longitude)
%%
% F_SUNRISE_SUNSET calculates the sunrise and sunset time of a geographic
% location based on the algorithm from Almanac for Computers, 1990. U.S.
% Government Printing Office.
%
% HOW: [day_night_map] = f_day_night_map(t_c,LATITUDE,LONGITUDE)
%
% Input: day_of_year    [1 x 1] number of days passes since Jan 0 of
%                               current year.
%
%        latitude       [n x m] latitude in degrees.
%
%        longitude      [n x m] longitude in degrees.
%
% Output: UT_sunrise    [n x m] UTC time of sunrise.
%
%         UT_sunset     [n x m] UTC time of sunset.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 23/01/2017
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 3
    error('Wrong number of input arguments')
end

if max(size(day_of_year)) ~= 1
    error('Select only one time epoch')
end

if max(size(latitude) ~= size(longitude) ) == 1
    error('Input locations do not have the same dimensions')
end


%Set correct limits for longitude
longitude(longitude > 180)            = longitude(longitude > 180) - 360;

%Calculate approximate time of phenomenon in days
t_sunrise                             = day_of_year + (6 - (longitude/15))/24;                                      %Days
t_sunset                              = day_of_year + (18 - (3*longitude/15))/24;                                   %Days

%Calculate Sun's zenith distance at rise or set
z                                     = 90 + (50/60);                                                               %Degrees

%Calculate Sun's mean anomaly
M_sunrise                             = 0.985600*t_sunrise - 3.289;                                                 %Degrees
M_sunset                              = 0.985600*t_sunset - 3.289;                                                  %Degrees

%Calculate Sun's true longitude
L_sunrise                             = M_sunrise + 1.916*sind(M_sunrise) + 0.020*sind(2*M_sunrise) + 282.634;      %Degrees
L_sunset                              = M_sunset + 1.916*sind(M_sunset) + 0.020*sind(2*M_sunset) + 282.634;         %Degrees

L_sunrise                             = mod(L_sunrise,360);                                                         %Degrees
L_sunset                              = mod(L_sunset,360);                                                          %Degrees

%Calculate Sun's right ascension
RA_sunrise                            = f_atan2d(0.91746*sind(L_sunrise),cosd(L_sunrise));                          %Degrees
RA_sunset                             = f_atan2d(0.91746*sind(L_sunset),cosd(L_sunset));                            %Degrees

%Calculate Sun's declination
sin_d_sunrise                         = 0.39782*sind(L_sunrise);
cos_d_sunrise                         = sqrt(1 - sin_d_sunrise.^2);

sin_d_sunset                          = 0.39782*sind(L_sunset);
cos_d_sunset                          = sqrt(1 - sin_d_sunset.^2);

%Calculate Sun's local hour angle
sin_H_sunrise                         = (cosd(z) - sin_d_sunrise.*sind(latitude))./(cos_d_sunrise.*cosd(latitude));
sin_H_sunset                          = (cosd(z) - sin_d_sunset.*sind(latitude))./(cos_d_sunset.*cosd(latitude));

index_sunrise_p                       = sin_H_sunrise > 1;
index_sunrise_m                       = sin_H_sunrise < -1;

index_sunset_p                        = sin_H_sunset > 1;
index_sunset_m                        = sin_H_sunset < -1;

sin_H_sunrise(abs(sin_H_sunrise) > 1) = NaN;
sin_H_sunset(abs(sin_H_sunset) > 1)   = NaN;

H_sunrise                             = 360 - acosd(sin_H_sunrise);                                                 %Degrees
H_sunset                              = acosd(sin_H_sunset);                                                        %Degrees

%Calculate local mean time of phenomenon
T_sunrise                             = H_sunrise/15 + RA_sunrise/15 - 0.065710*t_sunrise - 6.622;                  %Hours
T_sunset                              = H_sunset/15 + RA_sunset/15 - 0.065710*t_sunset - 6.622;                     %Hours

T_sunrise                             = mod(T_sunrise,24);                                                          %Hours
T_sunset                              = mod(T_sunset,24);                                                           %Hours

%Calculate universal time of phenomenon
UT_sunrise                            = day_of_year + (T_sunrise - longitude/15)/24;                                %Days
UT_sunset                             = day_of_year + (T_sunset - longitude/15)/24;                                 %Days

UT_sunrise(index_sunrise_p)           = +inf;
UT_sunrise(index_sunrise_m)           = -inf;

UT_sunset(index_sunset_p)             = +inf;
UT_sunset(index_sunset_m)             = -inf;

end
