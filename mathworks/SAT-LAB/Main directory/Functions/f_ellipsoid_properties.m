function [output] = f_ellipsoid_properties(name,type)
%%
% F_ELLIPSOID_PROPERTIES defines the values of the fundamental parameters
% of an ellipsoid.
%
% HOW: f = f_ellipsoid_properties(name)
%
% Input: name            String that specifies the reference ellipsoid
%                        name.
%
%        type            String that specifies the type of imported file.
%                        Should be 'xlsx' or 'csv'.
%
% Output: f             [1 x 1] element defined as ellipsoid.
%
% Dimitrios Piretzidis, Department of Rural and Surveying Engineering AUTh
% 2012
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin == 1
    
    type              = 'csv';
    
elseif nargin < 1 || nargin > 2
    
    error('Wrong number of input arguments')
    
end

output                = f_ellipsoid;
found_ell             = false;

if strcmp(type,'xlsx') == 1
    
    %Read the ellipsoid parameters using Ellipsoid.xlsx with Active X elements
    excel             = actxserver('Excel.Application');
    file              = excel.Workbooks.Open(which('Ellipsoid.xlsx'));
    sheet             = excel.Worksheets.get('Item', 'Ellipsoid');
    EPSG_data_range   = get(sheet,'Range', 'A2:H55');
    
    %Get the ellipsoid parameters
    EPSG_data         = EPSG_data_range.Value;
    
    %Close the excel file
    file.Close
    
    %Scan the data for the requested ellipsoid
    for i = 1:size(EPSG_data,1)
        
        if strcmp(name,EPSG_data{i,2}) == 1
            
            found_ell = true;
            data_row  = i;
            
        end
        
    end
    
elseif strcmp(type,'csv') == 1
    
    %Read the ellipsoid parameters using Ellipsoid.csv
    fid               = fopen('Ellipsoid.csv');
    
    %Get the ellipsoid parameters
    EPSG_header       = textscan(fid, '%s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');    
    EPSG_data         = textscan(fid, '%u %s %f %f %f %f %f %f', 54, 'Delimiter', ',');
    
    %Close the csv file
    fclose(fid);
    
    %Scan the data for the requested ellipsoid
    for i = 1:size(EPSG_data{1,2},1)
        
        if strcmp(name,EPSG_data{1,2}{i,1}) == 1
            
            found_ell = true;
            data_row  = i;
            
        end
        
    end
    
end

%Check if the requested ellipsoid has been found
if found_ell ~= true
    
    error('The requested ellipsoid has not been found. Check the ellipsoid name or choose different ellipsoid');
    
end

%Retrieve the ellipsoid parameters

if strcmp(type,'xlsx') == 1
    
    output.epsg_code  = EPSG_data{data_row,1};
    output.a          = EPSG_data{data_row,3};      %m
    output.b          = EPSG_data{data_row,4};      %m
    output.f          = EPSG_data{data_row,5};
    output.f_inv      = EPSG_data{data_row,6};
    output.e          = EPSG_data{data_row,7};
    output.e_prime    = EPSG_data{data_row,8};
    
elseif strcmp(type,'csv')==1
    
    output.epsg_code  = EPSG_data{1,1}(data_row,1);
    output.a          = EPSG_data{1,3}(data_row,1); %m
    output.b          = EPSG_data{1,4}(data_row,1); %m
    output.f          = EPSG_data{1,5}(data_row,1);
    output.f_inv      = EPSG_data{1,6}(data_row,1);
    output.e          = EPSG_data{1,7}(data_row,1);
    output.e_prime    = EPSG_data{1,8}(data_row,1);
    
end

output.gm             = 0.3986005E+15;              %m^3 s^-2
output.omega          = 7.292115E-05;               %rad/sec

end
