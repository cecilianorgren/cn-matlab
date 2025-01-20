function [data_ori,data_new] = f_read_tle(filename)
%%
% F_READ_TLE reads a NORAD two-line element (TLE) ephemeris file and
% extracts all the information for every satellite. For more information
% about the TLE format visit http://www.celestrak.com .
%
% HOW: [data_ori,data_new] = f_read_tle(filename)
%
% Input:  filename      [k x 1] string containing the TLE file name.
%
% Output: data_ori      [n x 1] structure containing all the TLE
%                               information in the original format.
%
%         data_new      [n x 1] structure containing all the TLE
%                               information in a processed format, ready to
%                               be used for orbital propagation.
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 16/12/2015
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 1
    error('Wrong number of input arguments')
end

%Load the file
file_ID                                             = fopen(filename,'r');

%Initialize the counter
i                                                   = 1;

%Read the file
while feof(file_ID) == 0
    
    current_line = fgetl(file_ID);
    
    if (strcmp(current_line(1,1:2),'1 ') == 0) && (strcmp(current_line(1,1:2),'2 ') == 0) && (isempty(current_line) == 0)
        
        data_ori(i,:).satellite_name                = current_line;         %Satellite Name
        data_new(i,:).satellite_name                = strtrim(current_line);
        
    elseif strcmp(current_line(1,1:2),'1 ') == 1
        
        %Get original data
        data_ori(i,:).line_element_1                = current_line(1,1);   %Line Number of Element Data
        data_ori(i,:).satellite_number_1            = current_line(3:7);   %Satellite NORAD Catalog Number
        data_ori(i,:).classification                = current_line(8);     %Classification (U=Unclassified)
        data_ori(i,:).launch_year                   = current_line(10:11); %International Designator (Last two digits of launch year)
        data_ori(i,:).launch_number                 = current_line(12:14); %International Designator (Launch number of the year)
        data_ori(i,:).launch_piece                  = current_line(15:17); %International Designator (Piece of the launch)
        data_ori(i,:).epoch_year                    = current_line(19:20); %Epoch Year (Last two digits of year)
        data_ori(i,:).epoch_day                     = current_line(21:32); %Epoch (Day of the year and fractional portion of the day)
        data_ori(i,:).first_derivative_mean_motion  = current_line(34:43); %First Time Derivative of the Mean Motion
        data_ori(i,:).second_derivative_mean_motion = current_line(45:52); %Second Time Derivative of Mean Motion (decimal point assumed)
        data_ori(i,:).BSTAR_drag_term               = current_line(54:61); %BSTAR drag term (decimal point assumed)
        data_ori(i,:).ephemeris_type                = current_line(63);    %Ephemeris type
        data_ori(i,:).element_number                = current_line(65:68); %Element number
        data_ori(i,:).checksum_1                    = current_line(69);    %Checksum (Modulo 10) (Letters, blanks, periods, plus signs = 0; minus signs = 1)
        
        %Process original data
        data_new(i,:).line_element_1                = str2double(current_line(1,1));
        data_new(i,:).satellite_number_1            = strtrim(current_line(3:7));
        data_new(i,:).classification                = strtrim(current_line(8));
        data_new(i,:).launch_year                   = str2double(current_line(10:11));
        
        if     data_new(i,:).launch_year > 56; data_new(i,:).launch_year = data_new(i,:).launch_year + 1900;
        elseif data_new(i,:).launch_year < 57; data_new(i,:).launch_year = data_new(i,:).launch_year + 2000;
        end
        
        data_new(i,:).launch_number                 = strtrim(current_line(12:14));
        data_new(i,:).launch_piece                  = strtrim(current_line(15:17));
        data_new(i,:).epoch_year                    = str2double(current_line(19:20));
        
        if     data_new(i,:).epoch_year > 56; data_new(i,:).epoch_year = data_new(i,:).epoch_year + 1900;
        elseif data_new(i,:).epoch_year < 57; data_new(i,:).epoch_year = data_new(i,:).epoch_year + 2000;
        end
        
        data_new(i,:).epoch_day                     = str2double(current_line(21:32));
        data_new(i,:).first_derivative_mean_motion  = 2*2*pi*str2double(current_line(34:43))/(86400^2);
        
        if     strcmp(current_line(45),'-') == 1                                     ; data_new(i,:).second_derivative_mean_motion = -str2double(['.',current_line(46:50)]);
        elseif strcmp(current_line(45),'+') == 1 || strcmp(current_line(45),' ') == 1; data_new(i,:).second_derivative_mean_motion =  str2double(['.',current_line(46:50)]);
        end
        
        data_new(i,:).second_derivative_mean_motion = data_new(i,:).second_derivative_mean_motion*10^str2double(current_line(51:52));
        data_new(i,:).second_derivative_mean_motion = 6*2*pi*data_new(i,:).second_derivative_mean_motion/(86400^3);
        
        if     strcmp(current_line(54),'-') == 1                                     ; data_new(i,:).BSTAR_drag_term = -str2double(['.',current_line(55:59)]);
        elseif strcmp(current_line(54),'+') == 1 || strcmp(current_line(54),' ') == 1; data_new(i,:).BSTAR_drag_term =  str2double(['.',current_line(55:59)]);
        end
        
        data_new(i,:).BSTAR_drag_term               = data_new(i,:).BSTAR_drag_term*10^(str2double(current_line(60:61)));
        data_new(i,:).ephemeris_type                = strtrim(current_line(63));
        data_new(i,:).element_number                = strtrim(current_line(65:68));
        data_new(i,:).checksum_1                    = strtrim(current_line(69));
        
    elseif strcmp(current_line(1,1:2),'2 ') == 1
        
        %Get original data
        data_ori(i,:).line_element_2                = current_line(1);     %Line Number of Element Data
        data_ori(i,:).satellite_number_2            = current_line(3:7);   %Satellite Number
        data_ori(i,:).inclination                   = current_line(9:16);  %Inclination [Degrees]
        data_ori(i,:).right_ascension_AN            = current_line(18:25); %Right Ascension of the Ascending Node [Degrees]
        data_ori(i,:).eccentricity                  = current_line(27:33); %Eccentricity (decimal point assumed)
        data_ori(i,:).argument_of_perigee           = current_line(35:42); %Argument of Perigee [Degrees]
        data_ori(i,:).mean_anomaly                  = current_line(44:51); %Mean Anomaly [Degrees]
        data_ori(i,:).mean_motion                   = current_line(53:63); %Mean Motion [Revs per day]
        data_ori(i,:).revolution_number             = current_line(64:68); %Revolution number at epoch [Revs]
        data_ori(i,:).checksum_2                    = current_line(69);    %Checksum (Modulo 10)
        
        %Process new data
        data_new(i,:).line_element_2                = strtrim(current_line(1));
        data_new(i,:).satellite_number_2            = strtrim(current_line(3:7));
        data_new(i,:).inclination                   = str2double(current_line(9:16));
        data_new(i,:).right_ascension_AN            = str2double(current_line(18:25));
        data_new(i,:).eccentricity                  = str2double(['.',current_line(27:33)]);
        data_new(i,:).argument_of_perigee           = str2double(current_line(35:42));
        data_new(i,:).mean_anomaly                  = str2double(current_line(44:51));
        data_new(i,:).mean_motion                   = 2*pi*str2double(current_line(53:63))/86400;
        data_new(i,:).revolution_number             = str2double(current_line(64:68));
        data_new(i,:).checksum_2                    = strtrim(current_line(69));
        
        %Update the counter
        i                                           = i + 1;
        
    end
    
end

%Close the file
fclose(file_ID);

end
