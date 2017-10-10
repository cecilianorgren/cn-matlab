dataset_id='C4_CP_EFW_L1_IB';
% On CFA
if 1
disp(char({' ','==== Testing CFA async ===='}))
URL_CFA=['http://cfadev.esac.esa.int/cfa/aio/async-product-action?USERNAME=avaivads&PASSWORD=!kjUY88lm&dataset_id=',dataset_id,'&START_DATE=2007-08-31T10:17:00.000000Z&END_DATE=2007-08-31T10:19:00.000000Z&DELIVERY_FORMAT=CDF&NON_BROWSER&NO_NOTIFY'];
disp(URL_CFA)
[fileName,st]=urlwrite(URL_CFA,fileName);
end
%
fid=fopen(fileName);
%
while 1    
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    disp(tline)
    if any(strfind(tline,'http:')) && any(strfind(tline,'gz')),       
        downloadfile = tline(strfind(tline,'http:'):strfind(tline,'gz')+1);
    end
end
%
tempDirectory=tempname;
mkdir(tempDirectory);
fileName=tempname;
gzFileName = [fileName '.gz'];
tic;[gzFileName,st]=urlwrite(downloadfile,gzFileName); cfat=toc;
gunzip(gzFileName);
%
fileNames=untar(fileName,[tempDirectory '/']);
for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

for iFile = 1:numel(fileNames)
    s = dir(fileNames{iFile});
    disp([num2str(s.bytes),' bytes downloaded in ',num2str(cfat),' seconds. Thats ',num2str(s.bytes*1e-6/cfat),' MB/second.'])
end