function pic2gif(fileDir,saveName,loops,delay)
% Adapted from code made by person below.
%
% This program creates a movie/slideshow from a set of images, and save it as an animated GIF file.
% Notice that the quality an image may decrease due to the GIF format.
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% September 2010 (C)


if strfind(saveName,'/')
    saveDir = '';
else
    saveDir = fileDir;
end
    
tmpFileNames = dir([fileDir '*.png']);
%tmpFileNames(1:2) = [];
nPics = numel(tmpFileNames);
fileName = {tmpFileNames.name}';
for k = 1:nPics
    %fileName{k} = tmpFileNames(k).name;        
    fileNum(k) = str2num(fileName{k}(regexp(fileName{k},'\d')));
end

[sortedFileNum,iSort]=sort(fileNum);
sortedFileName = fileName(iSort);

h = waitbar(0,['0% done'],'name','Progress') ;
for i=1:nPics  
    [M,c_map]=imread([fileDir,sortedFileName{i}],'png');
    if i==1
        imwrite(M,c_map,[saveDir,saveName],'gif','LoopCount',loops,'DelayTime',delay)
    else
        imwrite(M,c_map,[saveDir,saveName],'gif','WriteMode','append','DelayTime',delay)
    end    
    waitbar(i/length(sortedFileName),h,[num2str(round(100*i/length(sortedFileName))),'% done']) ;
end
close(h);
msgbox('Finished Successfully!')
disp(['Saved gif to ' saveDir ,saveName])
end
