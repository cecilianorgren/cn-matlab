function print(filename,varargin)
% CN.PRINT  Prints to "todays" directory.
%   CN.PRINT(filename,extrapath) checks if a directory under 
%   /Users/Cecilia/Research/<todays date here>
%   todays date is given by Matlab function datestr(now,'yyyy-mm-dd')
%   Saves as png, and adds extension .png to filename.

% If defined, construct additional path
frmt = 'png';
if nargin > 1
    if strcmp(varargin{1},'eps')
        frmt = 'epsc';
        extrapath='/';
    else
        extrapath = ['/',varargin{1},'/'];
    end
else
    extrapath='/';
end

% Make directory
directory = ['/amd/hem/export/home/cecilia/Research/ElectronHoleSimulation/Days/',datestr(now,'yyyy-mm-dd'),extrapath];

% If directory doesn'r already exist, create it
if ~exist(directory,'dir')
    eval(['mkdir ', directory])
end

% Add filename to directory path
path_and_file= [directory,filename,'.',frmt(1:3)];

% If a file with that name already exist, ask if you want to overwrite it
if exist(path_and_file,'file')
    question = 'File alreade exists, do you want to overwrite it? [1/0] >';
    print_flag = irf_ask(question,[],1);
else 
    print_flag=1;
end

% Prints in command window, what the action was
if print_flag
    set(gcf,'paperpositionmode','auto');
    eval(['print -d',frmt,' ',directory,filename,'.',frmt(1:3)])
    disp(['Printed ',path_and_file])
else
    disp('Did not print.')
end
    