function out = printpath()

extrapath = '/';

days = 1;
if days    
  [~,computername]=system('hostname');
  if strfind(computername,'ift0227887')
    directory = ['/Users/cno062/Research/Days/',datestr(now,'yyyy-mm-dd'),extrapath];
    
    directory_root = '/Users/cno062/Dropbox-IRFU/Cecilia Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
  elseif strfind(computername,'CeciliasMacBook')
    directory_root = '/Users/cecilia/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
    %directory = ['/Users/cecilia/Research/Days/',datestr(now,'yyyy-mm-dd'),extrapath];
  elseif strfind(computername,'Cecilias-MacBook')
    directory_root = '/Users/cecilia/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
    %directory = ['/Users/cecilia/Research/Days/',datestr(now,'yyyy-mm-dd'),extrapath];  
  elseif strfind(computername,'Cecilias-Mac-mini')
    directory_root = '/Users/cecilianorgren/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory_root = '/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
  elseif strfind(computername,'Mac')
    directory_root = '/Users/cecilianorgren/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory_root = '/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
  elseif strfind(computername,'IRFU018-cecilian.local')
    directory_root = '/Users/cecilianorgren/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory_root = '/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
  else % temporary
    directory_root = '/Users/cecilianorgren/IRFU\ Dropbox/Cecilia\ Norgren/Days/';
    directory_root = '/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/';
    directory = [directory_root, datestr(now,'yyyy-mm-dd'), extrapath];
  end 
  if ~exist(directory,'dir')
    %eval(['mkdir ', directory])
    mkdir(directory)
  end
end
out = directory;