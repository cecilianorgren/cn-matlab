function OK = setup_dir(Logical_file_id,root_save)

fileName = Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); 
numName = fileNameSplit{end-1};

dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
dirNameMatlab = sprintf('+event_%s%s%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
eventPath = ['/home/' localuser '/Research/Events/' dirName '/'];  
matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/' dirNameMatlab '/'];


mkdir(eventPath)