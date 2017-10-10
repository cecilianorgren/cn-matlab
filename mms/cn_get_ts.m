function [TS,dobj] = cn_get_ts(filePrefix,dataVar,tt)
% MMS.CN_GET_TS Load mms data file.
%   Simple loading routine that only checks that the file name should 
%   contain the given time interval. You need to verify that the file 
%   indeed contains the expected time interval manually.
%
%   [ts,dobj] = MMS.CN_GET_TS(filePrefix,fileVariable,time);
%   filePrefix - e.g. mms3_edp_brst_l2_scpot
%   fileVariable - e.g. mms3_edp_scpot
%   time - time in TSeries format
%
%   [ts,dobj] = MMS.CN_GET_TS(filePrefix,[],time);
%       Loads dataobject 'dobj' only.
%
%   MMS.CN_GET_TS(filePrefix,'list',time) - displays list of files in the
%       folder corresponding to the time given.
%

isDay = 0;
doLoad = 1;
doList = 0;

% Check if only a list is asked for.
if strcmp(dataVar,'list'); doLoad = 0; doList = 1; end

% Time representation
utc = tt.toUtc;
t.year  = str2double(utc(1:4));
t.month = str2double(utc(6:7));
t.day   = str2double(utc(9:10));
t.hour  = str2double(utc(12:13));
t.min   = str2double(utc(15:16));
t.sec   = str2double(utc(18:end-1));
tdateNum = datenum(t.year,t.month,t.day,t.hour,t.min,t.sec);
disp(['Time: ' utc])

% Directory information
dataDir = '/data/mms';
varDir = strjoin(strsplit(filePrefix,'_'),filesep);
dateDir = [utc(1:4) filesep utc(6:7) filesep utc(9:10)];
totDir = [dataDir filesep varDir filesep dateDir];

% List files in directory
listingD = dir([totDir filesep filePrefix '*.cdf']);
if isempty(listingD) % Files are in month folder
    totDir = totDir(1:end-3); isDay = 1;
    listingD = dir([totDir filesep filePrefix '*.cdf']);
    if isempty(listingD), % Also month folder is empty
        disp(['Could not find any ' filePrefix ' files during the same day!']); 
        TS = []; dobj = []; return;
    end
end
disp(['  File directory: ' totDir])
nFiles = numel(listingD);
if doList % Only list the files
    for ii = 1:nFiles, 
        disp(listingD(ii).name)        
    end
    TS = []; dobj = [];
    return;
end

% Find file that contains interval start time, t1
for ii = 1:nFiles,     
    splitName = strsplit(listingD(ii).name,'_');
    dateStr = splitName{end-1};
    if isDay, dateStr = [dateStr '000001']; end
    dateNums(ii,1) = datenum(str2double(dateStr(1:4)),str2double(dateStr(5:6)),str2double(dateStr(7:8)),str2double(dateStr(9:10)),str2double(dateStr(11:12)),str2double(dateStr(13:14))); 
end
% plot(1:nFiles,dateNums,'-*',[1 nFiles],tdateNum*[1 1])
fileId = find(tdateNum>dateNums,1,'last');
if isempty(fileId); 
    fileId = 1; 
    allId = fileId;
else
    lastId = fileId;
    allId = find(dateNums==dateNums(lastId));
end
nId = numel(allId);


if nId > 1            
    %specId = input(['Choose file #(' num2str(1) '-' num2str(nId) '): ']);
    %fileId = allId(specId);
    fileId = allId(end); % just load latest version
    %disp('Loading latest version.')    
end
disp(['  File: ' listingD(fileId).name ' (' num2str(listingD(fileId).bytes) ' bytes)'])

%% Load the data
dobj=dataobj([totDir filesep listingD(fileId).name]);
if isempty(dataVar)
    TS = [];
else    
    if ischar(dataVar)
        dataVar = {dataVar};
    end
    for ii = 1:numel(dataVar)
        try   
            ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));            
        catch    
            ts = get_variable(dobj,dataVar{ii});            
        end    
        TS(ii) = ts;
    end
end