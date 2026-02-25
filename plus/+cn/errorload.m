loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/';
allMats = ls([loadPath 'Spis/*.mat']);
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
if plotAll; plotThis = 1:nFiles; else plotThis = nFiles; end

isError = 1;
while isError
    try
        ExB.plot1;        
    catch errorMsg                
        isError = 0;
    end
end