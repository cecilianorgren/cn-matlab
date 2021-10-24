function print(filename,varargin)
% CN.PRINT  Prints to "todays" directory.
%   CN.PRINT(filename,options) checks if a directory under 
%   /Users/Cecilia/Research/<todays date here>
%   todays date is given by Matlab function datestr(now,'yyyy-mm-dd')
%   Saves as png, and adds extension .png to filename.
%   
%   options

% Defaults

filepath = '/';
frmt = 'png';
fileending = 'png';
print_flag = 1;
have_options = 0;
extrapath = '/';
doIncrement = 0;
days = 1;
printInLocal = 0;
renderer = '-painters';
printInPlotTime  = 0;
isharddirectory = 0;

% Collect options
if nargin > 1; have_options = 1; end  
while have_options
    %if ischar(varargin{1});
    %    if findstr(varargin{1})
    %        extrapath = ['/',varargin{1},'/'];
    %        varargin(1) = [];
    %    end
    %end
    switch lower(varargin{1})
      case 'eps'
          frmt = 'epsc';
          fileending = 'eps';
          varargin(1) = [];
      case 'jpeg'
          frmt = 'jpeg';
          fileending = 'jpg';
          varargin(1) = [];
      case 'pdf'
          frmt = 'pdf';
          fileending = 'pdf';
          varargin(1) = [];
      case 'incr'
          doIncrement = 1;    
      case 'hard_path'
        isharddirectory = 1;        
        harddirectory = varargin{2};
        doIncrement = 1;               
        varargin(1:2) = [];  
      case 'path'
          days = 0;
          directory = varargin{2};
          doIncrement = 1;               
          varargin(1:2) = [];    
      case 'here'
          days = 0;
          directory = pwd; directory = [pwd '/'];
          varargin(1) = []; 
      case 'thor'
        isharddirectory = 1;
        harddirectory = '/Users/Cecilia/Research/THOR/';
        varargin(1) = [];
      case 'time'
        printInPlotTime = 1;
        time = varargin{2};
        varargin(1) = [];
        varargin(1) = [];
      case 'painters'
        renderer = '-painters';
        varargin(1) = [];
      case 'opengl'
        renderer = '-opengl';
        varargin(1) = [];
      otherwise
        disp(sprintf(' -- Can not recognize input %s',varargin{1}))
        varargin{1} = [];
    end
    if isempty(varargin); break; end
end

% If defined, construct additional path
% if nargin > 1    
%     if strcmp(varargin{1},'eps')
%         frmt = 'epsc';
%         extrapath='/';
%     else
%         extrapath = ['/',varargin{1},'/'];
%     end
% else
%     extrapath='/';
% end

% Make directory. If directory doesn'r already exist, create it
if isharddirectory
  directory = harddirectory;
elseif printInLocal
  directory = '';
elseif printInPlotTime
  timeutc = time.utc;
  [~,computername]=system('hostname');
  if strfind(computername,'ift0227887')
    directory = ['/Users/cno062/Research/Events/',timeutc(1:10),extrapath];
  else
    directory = ['/Users/Cecilia/Research/Events/',timeutc(1:10),extrapath];
  end  
  
  if ~exist(directory,'dir')
    eval(['mkdir ', directory])
  end
else
  if days    
    [~,computername]=system('hostname');
    if strfind(computername,'ift0227887')
      directory = ['/Users/cno062/Research/Days/',datestr(now,'yyyy-mm-dd'),extrapath];
    else
      directory = ['/Users/Cecilia/Research/Days/',datestr(now,'yyyy-mm-dd'),extrapath];
    end 
    if ~exist(directory,'dir')
      eval(['mkdir ', directory])
    end
  end
end
% Add filename to directory path
%path_and_file= [directory,filename,'.',frmt(1:3)];
path_and_file= [directory,filename,'.',fileending];


% Check if file already exists-
% If a file with that name already exist, 
% ask if you want to overwrite it, or make increment
%if doIncrement
%    
%else
if exist(path_and_file,'file')
    question = 'File already exists, do you want to overwrite it? [1/0] >';
    print_flag = irf_ask(question,[],1);
%else 
    
end

% Prints in command window, what the action was
if print_flag
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');
    %set(gcf,'color','none');
    %eval(['print -d',frmt,' ',directory,filename,'.',frmt(1:3)])
    %print(['-d',frmt],'-opengl','-r300',[directory,filename,'.',frmt(1:3)]);
    %print(['-d',frmt],'-painters','-r300',[directory,filename,'.',frmt(1:3)]);
    %print(['-d',frmt],'-painters','-r600',[directory,filename,'.',frmt(1:3)]);
    %print(['-d',frmt],'-opengl','-r300',[directory,filename,'.',frmt(1:3)]);
    print(['-d',frmt],renderer,'-r300',[directory,filename,'.',fileending]);
    disp(['Printed ',path_and_file])
else
    disp('Did not print.')
end
    