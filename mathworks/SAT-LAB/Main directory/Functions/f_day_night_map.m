function [day_night_map] = f_day_night_map(t_c,LATITUDE,LONGITUDE)
%%
% F_DAY_NIGHT_MAP returns the day/night map of the input geographic
% locations.
%
% HOW: [day_night_map] = f_day_night_map(t_c,LATITUDE,LONGITUDE)
%
% Input: t_c            [1 x 1] UTC time in the form of day of year.
%
%        LATITUDE       [n x m] latitude in degrees.
%
%        LONGITUDE      [n x m] longitude in degrees.
%
% Output: day_night_map [n x m] matrix with NaNs in the locations where is
%                               still day and -100000 in the locations 
%                               where is still night.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 23/01/2017
%
% uses m-files: f_sunrise_sunset

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 3
    error('Wrong number of input arguments')
end

if max(size(t_c)) ~= 1
    error('Select only one time epoch')
end

if max(size(LATITUDE) ~= size(LONGITUDE) ) == 1
    error('Input locations do not have the same dimensions')
end

%Get day of year
t_0                                      = floor(t_c);

%Calculate sunrise-sunset time
[R,S]                                    = f_sunrise_sunset(t_0,LATITUDE,LONGITUDE);

%Create map layer
day_night_map                            =-100000*ones(size(LONGITUDE));

%Set day points
day_night_map(t_c > R & t_c < S)         = 100000;
day_night_map(t_c - 1 > R | t_c + 1 < S) = 100000;
day_night_map(R == -inf)                 = 100000;
day_night_map(R == +inf)                 =-100000;

%Fix some inconsistences
for i = 2:size(day_night_map,1) - 1
    
    for j = 2:size(day_night_map,2) - 1
        
        neighbors                        = [day_night_map(i-1,j-1),day_night_map(i-1,j),day_night_map(i-1,j+1),...
                                            day_night_map(i,j-1)                       ,day_night_map(i,j+1)  ,...
                                            day_night_map(i+1,j-1),day_night_map(i+1,j),day_night_map(i+1,j+1)];
        
        neighbors_sum                    = sum(neighbors);
        
        if sign(neighbors_sum) ~= sign(day_night_map(i,j))
            
            day_night_map(i,j)           = -day_night_map(i,j);
            
        end
        
    end
    
end

day_night_map(day_night_map == 100000)   = NaN;

end
