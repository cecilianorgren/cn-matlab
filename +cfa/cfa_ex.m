%tmpDir = tempname;
%mkdir(tmpDir);
%cd(tmpDir);

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

%temp_file='delme';
%disp('== Download:')
%urlDownloadLine='http://cfadev.esac.esa.int/cfa/aio/product-action?&NON_BROWSER&DATASET_ID=C3_CP_ASP_IONC&START_DATE=2004-06-18T11:35:00Z&END_DATE=2004-06-19T18:35:00Z';
%[downloadMsg,downloadStatus]=urlwrite(urlDownloadLine,temp_file);
%disp(urlDownloadLine)
%disp(downloadMsg(1:end))

disp('== Logout:')
urlLogoutLine='http://cfadev.esac.esa.int/cfa/aio/login-action?SUBMIT=LOGOUT&NON_BROWSER';
[logoutMsg,logoutStatus]=urlread(urlLogoutLine,'timeout',10);
disp(urlLogoutLine)
disp(logoutMsg(1:end-1))

%% Arnauds ex
URL = 'http://cfadev.esac.esa.int/cfa/aio/product-action';
Data = urlread(URL, 'Authentication', 'Basic','Get', { 'Username', 'avaivads', 'password', '!kjUY88lm', 'DATASET_ID', 'C3_CP_ASP_IONC', 'START_DATE', '2004-06-18T11:35:00Z', 'END_DATE', '2004-06-19T18:35:00Z', 'NON_BROWSER', '1'});