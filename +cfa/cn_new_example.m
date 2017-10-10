%tmpDir = tempname;
%mkdir(tmpDir);
cd('/Users/Cecilia/CFA');

if 1, % Example 1 % 2 hours of data
    disp('==== Example 1: C3_CP_RAP_I3D_H ====')
    disp('==== CFA ====')
    tic;cfa.download('2004-06-18T11:35:00Z/2004-06-19T18:35:00Z','C3_CP_ASP_IONC','nowildcard'); toc 
    %delete(temp_file); 
    disp('==== CAA ====')
    % ********** CAA **************
    %tic;caa_download('2001-02-01T00:00:00.000Z/2001-02-01T02:00:00.000Z','C3_CP_RAP_I3D_H','nowildcard'); toc  
end
%rmdir(tmpDir,'s');
%%
tmpDir = tempname;
mkdir(tmpDir);
cd(tmpDir);

disp(char({' ','==== Testing CFA ===='}))
disp('== Login:')
urlLoginLine='http://cfadev.esac.esa.int/cfa/aio/login-action?USERNAME=avaivads&PASSWORD=!kjUY88lm&SUBMIT=LOGIN&NON_BROWSER';
[loginMsg,loginStatus]=urlread(urlLoginLine,'timeout',10);
disp(urlLoginLine)
disp(loginMsg(1:end-1))

disp('== Download:')
urlDownloadLine='http://cfadev.esac.esa.int/cfa/aio/product-action?&NON_BROWSER&DATASET_ID=C3_CP_ASP_IONC&START_DATE=2004-06-18T11:35:00Z&END_DATE=2004-06-19T18:35:00Z';
[downloadMsg,downloadStatus]=urlread(urlDownloadLine);
disp(urlDownloadLine)
disp(downloadMsg(1:end-1))

temp_file='delme';
disp('== Download:')
urlDownloadLine='http://cfadev.esac.esa.int/cfa/aio/product-action?&NON_BROWSER&DATASET_ID=C3_CP_ASP_IONC&START_DATE=2004-06-18T11:35:00Z&END_DATE=2004-06-19T18:35:00Z';
[downloadMsg,downloadStatus]=urlwrite(urlDownloadLine,temp_file);
disp(urlDownloadLine)
disp(downloadMsg(1:end))

disp('== Logout:')
urlLogoutLine='http://cfadev.esac.esa.int/cfa/aio/login-action?SUBMIT=LOGOUT&NON_BROWSER';
[logoutMsg,logoutStatus]=urlread(urlLogoutLine,'timeout',10);
disp(urlLogoutLine)
disp(logoutMsg(1:end-1))
%disp('== Download:')
%urlDownloadLine='http://cfadev.esac.esa.int/cfa/aio/product-action?&DATASET_ID=C3_CP_ASP_IONC&START_DATE=2004-06-18T11:35:00Z&END_DATE=2004-06-19T18:35:00Z';
%[downloadMsg,downloadStatus]=urlread(urlDownloadLine);
%disp(urlDownloadLine)
%disp(downloadMsg(1:end-1))

%%
temp_file='delme';
disp('== Download:')
urlDownloadLine='http://cfadev.esac.esa.int/cfa/aio/product-action?USERNAME=avaivads&PASSWORD=!kjUY88lm&DATASET_ID=C3_CP_ASP_IONC&START_DATE=2004-06-18T11:35:00Z&END_DATE=2004-06-19T18:35:00Z&NON_BROWSER';
[downloadMsg,downloadStatus]=urlwrite(urlDownloadLine,temp_file);
disp(urlDownloadLine)
disp(downloadMsg(1:end))
