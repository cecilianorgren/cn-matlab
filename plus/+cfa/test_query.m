% Test query requests
tint = toepoch([2007 08 31 10 17 00;2007 08 31 10 19 00])';
isoT1 = cfa.time(tint(1)); isoT2 = cfa.time(tint(2));
%%
%isoT1 = '2004-12-06T23:50:00Z';
%isoT2 = '2004-12-07T01:00:00Z';
cfaServer = 'http://cfadev.esac.esa.int/';
cfaAction = 'cfa/aio/metadata-action?';
cfaIdentity='&USERNAME=avaivads&PASSWORD=!kjUY88lm';
datasetId = 'C%_CP_EFW_L3_E3D_INERT';
datasetId=urlencode(datasetId);
%datasetId = 'C%2z_EDI_%25_CDF';
selectedFields = 'SELECTED_FIELDS=DATASET_INVENTORY&RESOURCE_CLASS=DATASET_INVENTORY';
nonBrowser = '&NON_BROWSER';
wildcard = '%25'; % %-sign
space = '%20'; % space-sign
and = [space 'AND' space];
fileName = [tempname '.csv'];
%and = [' AND '];
%formatVOTable = 'VOTABLE';
%formatJSON = 'JSON';
%formatCSV = 'CSV';
%format={'VOTABLE','JSON','CSV'};
%queryFormat = ['&RETURN_TYPE=' format];
queryFormat = {'&RETURN_TYPE=VOTABLE','&RETURN_TYPE=JSON','&RETURN_TYPE=CSV'};
queryData = ['&QUERY=DATASET.DATASET_ID' space 'like' space '''' datasetId ''''];
queryTime = [and 'DATASET_INVENTORY.START_TIME' space '<=' space '''' isoT2 '''',...
             and 'DATASET_INVENTORY.END_TIME' space '>=' space '''' isoT1 ''''];
queryLine = [cfaServer cfaAction selectedFields,...
             queryData queryTime queryFormat{3} '&NO_NOTIFY'];
disp(queryLine)

caalog=urlread(queryLine)
disp(caalog)


caafile=urlwrite(queryLine,fileName)
csvread(caafile)

%xslt(caafile)
%xDoc = xmlread(caafile);
%% Tidy up caalog and disp
caa = textscan(caalog,'%s%s%s%s%s','delimiter',',');
nDataset = numel(caa{1})-1;
nCharDataset = 40;
nCharData = 24;
nCharInstances = 12;
nCharVersion = 4;
caaout = '"DATASET.DATASET_ID","DATASET_INVENTORY.START_TIME","DATASET_INVENTORY.END_TIME","DATASET_INVENTORY.NUM_INSTANCES","DATASET_INVENTORY.VERSION"';
for ii = 1:nDataset
    %caaout = 1
end
    


%% try urlread with example from manual
%manualLine='http://cfadev.esac.esa.int/cfa/aio/metadata-action?SELECTED_FIELDS=DATASET_INVENTORY&RESOURCE_CLASS=DATASET_INVENTORY&QUERY=DATASET.DATASET_ID%20like%20''C%25_EDI_%25_CDF''%20AND%20DATASET_INVENTORY.START_TIME%20<=%20''2004-12-07T01:00:00Z''%20AND%20DATASET_INVENTORY.END_TIME%20>=%20''2004-12-06T23:50:00Z''&NONBROWSER';
manualLine='http://cfadev.esac.esa.int/cfa/aio/metadata-action?SELECTED_FIELDS=DATASET_INVENTORY&RESOURCE_CLASS=DATASET_INVENTORY&QUERY=DATASET.DATASET_ID%20like%20''C%25_EDI_%25_CDF''%20AND%20DATASET_INVENTORY.START_TIME%20<=%20''2004-12-07T01:00:00Z''%20AND%20DATASET_INVENTORY.END_TIME%20>=%20''2004-12-06T23:50:00Z''';
disp(manualLine)
disp(queryLine)
disp(['Are they identical? ' num2str(strcmp(manualLine,queryLine))]);
% no they weren't :(, but now they are! ^
urlread(manualLine)
urlread(queryLine)