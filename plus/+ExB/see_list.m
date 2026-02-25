function out = see_list(varargin)

loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/';
allMats = ls([loadPath 'Spis/*.mat']);

if nargout == 0 % no output, just display list   
    disp(allMats)
else % some sort of flag, make cell list
    indexDot = strfind(allMats,'.');
    indexStart = [1 indexDot(1:(end-1))+5];
    indexEnd = indexDot+3;
    indexFileSep = strfind(allMats,'/');
    nFiles = numel(indexDot);
    vecToLoad = cell(1,nFiles);
    vecFileName = cell(1,nFiles);
    for kk = 1:nFiles, vecToLoad{kk} = allMats(indexStart(kk):indexEnd(kk)); end
    for pp = 1:nFiles
        son = find(indexFileSep<indexEnd(pp),1,'last'); 
        vecFileName{pp} = allMats(indexFileSep(son)+1:indexEnd(pp)); 
    end
    % check for flag
    if nargin ~= 0
    flag = varargin{1};    
    switch flag
        case 'last'
            if nargin == 1
                out = vecFileName(end);            
            else
                nOut = varargin{2};
                out = vecFileName((end-nOut+1):end);
            end
        case 'first'
            if nargin == 1
                out = vecFileName(1);            
            else
                nOut = varargin{2};
                out = vecFileName(1:nOut);
            end
        case 'all'
            out = vecFileName(:);
    end            
    else out = vecFileName(:);
    end
end
            