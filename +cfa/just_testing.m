dataset_id='C*_CP_WHI_ELECTRON_DENSITY';
% On CFA
URL_CFA=['http://cfadev.esac.esa.int/cfa/aio/product-action?USERNAME=avaivads&PASSWORD=!kjUY88lm&dataset_id=',dataset_id,'&START_DATE=2010-01-01T00:00:00Z&END_DATE=2010-01-01T23:59:59Z&NON_BROWSER&NO_NOTIFY'];
fileName=tempname;
gzFileName = [fileName '.tar.gz'];
tarFileName = [fileName '.tar'];
[gzFileName,st]=urlwrite(URL_CFA,gzFileName); 
if 1
    tic 
    gunzip(gzFileName)
    fileNames=untar(fileName);
    toc
else
    tic
    fileNames=untar(gzFileName);
    toc
end
for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end



%% On CAA
disp(char({' ','==== Testing CAA ===='}))
URL_CAA=['http://caa.estec.esa.int/caa_query?uname=Norgren01&pwd=sally&dataset_id=',dataset_id,'&time_range=2010-01-01T00:00:00Z/2010-01-01T23:59:59Z&nonotify=1&refdoc=0'];
disp(URL_CAA)
fileName=tempname;
zFileName = [fileName '.zip'];
[zFileName,st]=urlwrite(URL_CAA,zFileName); 
fileNames=unzip(zFileName);
for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

%%
fileName=tempname;
zFileName = [fileName '.zip'];
gzFileName = [fileName '.gz'];

[dl_zFileName,stCAA]=urlwrite(URL_CAA,fileName); 
[dl_gzFileName,stCFA]=urlwrite(URL_CFA,gzFileName); 
[dl_fileName,stCFAb]=urlwrite(URL_CFA,fileName); 
%%
disp(dl_zFileName)
disp(dl_gzFileName)
disp(dl_fileName)
%%
caa=unzip(dl_zFileName);
cfa_1=untar(dl_gzFileName);
cfa_2=untar(dl_fileName);
%%
disp(caa{1})
disp(cfa_1{1})
disp(cfa_2{1})


%fileNames=unzip(zFileName);
