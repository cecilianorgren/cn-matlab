function [t s] = dl(dataset_id)

%% On CFA
try    
    identityCFA = 'USERNAME=avaivads&PASSWORD=!kjUY88lm';
    urlNonotify = '&NO_NOTIFY';
    disp(char({' ','==== Testing CFA ===='}))
    URL_CFA=['http://cfadev.esac.esa.int/cfa/aio/product-action?',identityCFA,'&dataset_id=',dataset_id,'&START_DATE=2010-01-01T00:00:00Z&END_DATE=2010-01-01T23:59:59Z&NON_BROWSER',urlNonotify];
    disp(URL_CFA)
    fileName=tempname;
    gzFileName = [fileName '.gz'];
    tic;[gzFileName,st]=urlwrite(URL_CFA,gzFileName); tCFA=toc;
    gunzip(gzFileName);
    fileNames=untar(fileName);
    for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

    for iFile = 1:numel(fileNames)
        sCFA = dir(fileNames{iFile});
        disp([num2str(sCFA.bytes),' bytes downloaded in ',num2str(tCFA),' seconds. Thats ',num2str(sCFA.bytes*1e-6/tCFA),' MB/second.'])
    end
    
    t.CFA=tCFA;
    s.CFA=sCFA.bytes;
catch 
    disp('Error downloading CFA data.')    
    t.CFA=0;
    s.CFA=0;
end

%% On CAA
try
    identityCAA='uname=Norgren01&pwd=sally';
    urlNonotify = '&nonotify=1';
    disp(char({' ','==== Testing CAA ===='}))
    URL_CAA=['http://caa.estec.esa.int/caa_query?',identityCAA,'&dataset_id=',dataset_id,'&time_range=2010-01-01T00:00:00Z/2010-01-01T23:59:59Z&refdoc=0',urlNonotify];
    disp(URL_CAA)
    fileName=tempname;
    zFileName = [fileName '.zip'];
    tic;[zFileName,st]=urlwrite(URL_CAA,zFileName); tCAA=toc;
    fileNames=unzip(zFileName);
    for iFile = 1:numel(fileNames), disp(fileNames{iFile}), end

    for iFile = 1:numel(fileNames)
        sCAA = dir(fileNames{iFile});
        disp([num2str(sCAA.bytes),' bytes downloaded in ',num2str(tCAA),' seconds. Thats ',num2str(sCAA.bytes*1e-6/tCAA),' MB/second.'])
    end

    t.CAA=tCAA;    
    s.CAA=sCAA.bytes;
    
catch 
    disp('Error downloading CAA data.')    
    t.CAA=0;
    s.CAA=0;
end

%% Compare times
disp(char({' ',['===> t_CAA/t_CFA=',num2str(t.CAA/t.CFA)]}));