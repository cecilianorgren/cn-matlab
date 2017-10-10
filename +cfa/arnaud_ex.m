dataset_id='C4_CP_EFW_L1_IB';
%% On CFA
disp(char({' ','==== Testing CFA ===='}))
URL_CFA=['http://cfadev.esac.esa.int/cfa/aio/async-product-action?USERNAME=avaivads&PASSWORD=!kjUY88lm&dataset_id=',dataset_id,'&START_DATE=2007-08-31T10:17:00.000000Z&END_DATE=2007-08-31T10:19:00.000000Z&DELIVERY_FORMAT=CDF&NON_BROWSER&NO_NOTIFY'];
disp(URL_CFA)
tempDirectory=tempname;
mkdir(tempDirectory);
fileName=tempname;
gzFileName = [fileName '.gz'];
tic;[gzFileName,st]=urlwrite(URL_CFA,gzFileName); cfat=toc;
gunzip(gzFileName);
fileNames=untar(fileName,[tempDirectory '/']);
for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

for iFile = 1:numel(fileNames)
    s = dir(fileNames{iFile});
    disp([num2str(s.bytes),' bytes downloaded in ',num2str(cfat),' seconds. Thats ',num2str(s.bytes*1e-6/cfat),' MB/second.'])
end

%% On CAA
disp(char({' ','==== Testing CAA ===='}))
URL_CAA=['http://caa.estec.esa.int/caa_query?uname=vaivads&pwd=caa&dataset_id=',dataset_id,'&time_range=2010-01-01T00:00:00Z/2010-01-01T23:59:59Z&refdoc=0'];
disp(URL_CAA)
fileName=tempname;
zFileName = [fileName '.zip'];
tic;[zFileName,st]=urlwrite(URL_CAA,zFileName); caat=toc;
fileNames=unzip(zFileName);
for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

for iFile = 1:numel(fileNames)
    s = dir(fileNames{iFile});
    disp([num2str(s.bytes),' bytes downloaded in ',num2str(caat),' seconds. Thats ',num2str(s.bytes*1e-6/caat),' MB/second.'])
end

%% Compare times

disp({' ',['===> t_CAA/t_CFA=',num2str(caat/cfat)]})