% Merge different run who has the same parameters.
toMerge = {'Spis/20140409T222215_Tao2011_3051082.mat',...
           'Spis/20140409T222215_Tao2011_3051082.mat'}
loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/';       

meanVxyzTmp = zeros(size(meanVxyzTmp));
meanVXxyzTmp = zeros(size(meanVXxyzTmp));
meanVYxyzTmp = zeros(size(meanVYxyzTmp));
meanVZxyzTmp = zeros(size(meanVZxyzTmp));

sumMxyzTmp = zeros(size(sumMxyzTmp));
sumMVxyxTmp = zeros(size(sumMVxyxTmp));
sumMVXxyxTmp = zeros(size(sumMVXxyxTmp));
sumMVYxyxTmp = zeros(size(sumMVYxyxTmp));
sumMVZxyxTmp = zeros(size(sumMVZxyxTmp));

sumMVrzTmp = zeros(size());
sumMrzTmp = zeros(size());
       

% First load everything for one run
load([loadPath toMerge{1}]);

% Then load the other runs and add to appropriate array/matrix
for ii=2:numel(toMerge)
    sumMxyx = sumMxyx + load([loadPath toMerge{ii}],'sumMxyz');
    sumVxyx = sumVxyx + load([loadPath toMerge{ii}],'sumVxyz');
    sumVXxyx = sumVXxyx + load([loadPath toMerge{ii}],'sumVXxyz');
    sumVYxyx = sumVYxyx + load([loadPath toMerge{ii}],'sumVYxyz');
    sumVZxyx = sumVZxyx + load([loadPath toMerge{ii}],'sumVZxyz');
    
    variables{ii} = load([loadPath toMerge{ii}],...
        'sumMxyz','sumMVxyx','sumMVXxyx','sumMVYxyx','sumMVZxyx','x0','y0','z0','xc','yc','',''
    