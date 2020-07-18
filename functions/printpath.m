function out = printpath()

extrapath = '/';

days = 1;
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
out = directory;