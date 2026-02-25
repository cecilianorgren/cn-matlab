% Testing CFA vs CAA.
doPrint = 1;
startTimeStr = datestr(now,'yyyymmddTHHMMSS');
for doStream = [0 1];
if doStream, streamStr = 'stream'; else streamStr = 'nostream'; end
fileDir = '/Users/Cecilia/CFA/';
fileName = ['test_caa_download_',startTimeStr,'_',streamStr '.txt'];
filepath = [fileDir fileName];
fid = fopen(filepath,'w');
header = {'Dataset                                 Start time                      End time                        CAA Bytes       CSA Bytes       CAA time        CSA time        t_CSA/t_CAA',...
          '***********                             ***********                     ***********                     ***********     ***********     ***********     ***********     ***********'};
fprintf(fid,[header{1},'\n',header{2},'\n']);
fclose(fid);
%%

fid = fopen(filepath,'a');

% Datasets to download
dataset={'C3_CP_EFW_L3_E3D_INERT',...
         'C3_CP_PEA_3DXPH_cnts',...
         'C3_CP_RAP_PAD_L3DD',...
         'C3_CP_WHI_ELECTRON_DENSITY',...
         'C3_CP_STA_DWF_NBR',...
         'C2_CP_RAP_EPITCH',...
         'C3_CP_RAP_DE',...
         'C3_CP_PEA_PADMARH_DEFlux',...
         'C3_CP_FGM_FULL',...
         'C3_CP_EFW_L2_V3D_GSE_EX',...                          
         'C3_CP_CIS-CODIF_HS_H1_MOMENTS',...
         'C4_CP_CIS-CODIF_RPA_HE1_PSD',...
         'C2_CP_PEA_PITCH_FULL_PSD',...
         'C3_CP_WHI_ELECTRON_GYROFREQUENCY',...         
         'C3_CP_STA_PSD',...
         'C3_CP_ASP_IONC',...
         'C1_CP_DWP_CORR_ST',...
         'C1_CP_EDI_EGD',...
         'C4_CG_WBD_GIFPLOT',...
         }; 
numberOfDatasets = numel(dataset);

% Iterate through datasets and writes results to file
% Dataset is 40 char
% Start time / End time is 32 char
% Bytes is 16 char
% CFA/CAA download time is 24
% time ratio is 16 char
procedure = {'simple','caa_download'};
t = cell(numberOfDatasets,1);
s = cell(numberOfDatasets,1);
for ii=1:numberOfDatasets
    % Time interval, new for each dataset
    t1 = [2005 01 29 floor(rand(1,1)*24) floor(rand(1,1)*60) 00]; % sometime on Jan 29th
    tT = 2+floor(rand(1,1)*24); % duration, somewhere between 2 and 24 hours
    tint = toepoch(t1)+[0 tT*60*60];
    
    switch procedure{2}
        case 'simple'
            [t{ii} s{ii}] = cfa.dl(dataset{ii});
        case 'caa_download'
            try % downloading from CSA
                tic; 
                if doStream; caa_download(tint,dataset{ii},'downloadDirectory=test_csa','stream','nowildcard'); 
                else         caa_download(tint,dataset{ii},'downloadDirectory=test_csa','nowildcard'); 
                end
                t{ii}.CSA = toc; 
                dirs = dir(['test_csa/' dataset{ii}]); 
                s{ii}.CSA = dirs(3).bytes;                
            catch 
                disp('Error downloading CSA data.')    
                t{ii}.CSA=0;
                s{ii}.CSA=0;
            end
            try  % downloading from CAA
                tic; 
                if doStream caa_download(tint,dataset{ii},'caa','downloadDirectory=test_caa','stream','nowildcard'); 
                else        caa_download(tint,dataset{ii},'caa','downloadDirectory=test_caa','nowildcard'); 
                end
                t{ii}.CAA = toc;
                dirs = dir(['test_caa/' dataset{ii}]); 
                s{ii}.CAA = dirs(3).bytes;
            catch 
                disp('Error downloading CAA data.')    
                t{ii}.CAA=0;
                s{ii}.CAA=0;
            end
    end     
       
    % print to text file
    eval(['str{1} = sprintf(''%s %',num2str(40-numel(dataset{ii})-1),'s'',dataset{ii},'' '');'])   
    isoT1 = irf_time(tint(1),'epoch2iso');
    isoT2 = irf_time(tint(2),'epoch2iso'); 
    eval(['str{2} = sprintf(''%s %',num2str(32-numel(isoT1)-1),'s'',isoT1,'' '');'])
    eval(['str{3} = sprintf(''%s %',num2str(32-numel(isoT2)-1),'s'',isoT2,'' '');'])
    s4=sprintf('%.0f',s{ii}.CAA);
    eval(['str{4} = sprintf(''%s % ',num2str(16-numel(s4)-1),'s'',s4,'' '');'])
    s5=sprintf('%.0f',s{ii}.CSA);
    eval(['str{5} = sprintf(''%s % ',num2str(16-numel(s5)-1),'s'',s5,'' '');'])
    s6=sprintf('%f',t{ii}.CAA);
    eval(['str{6} = sprintf(''%s % ',num2str(16-numel(s6)-1),'s'',s6,'' '');'])
    s7=sprintf('%f',t{ii}.CSA);
    eval(['str{7} = sprintf(''%s % ',num2str(16-numel(s7)-1),'s'',s7,'' '');'])
    s8=sprintf('%f',t{ii}.CSA/t{ii}.CAA); 
    eval(['str{8} = sprintf(''%s % ',num2str(16-numel(s8)-1),'s'',s8,'' '');'])
    str{9} = '\n'; % new line
    
    fprintf(fid,[str{1},str{2},str{3},str{4},str{5},str{6},str{7},str{8},str{9}]);  
end
fclose('all');

%% Read data from text file and present in plot time ratio vs time
fid = fopen(filepath,'r');
textFormat = '%s%s%s%f%f%f%f%f';
data = textscan(fid,textFormat,'headerlines',2);
fclose(fid);

cellt1 = data{2}; 
cellt2 = data{3}; 
for oo = 1:numel(cellt1); 
    t1(oo,1) = iso2epoch(cellt1{oo}); 
    t2(oo,1) = iso2epoch(cellt2{oo}); 
end
% t1 = deal(data{2});
% t2 = deal(data{3});
%
%
bytesCAA = data{4};
bytesCSA = data{5};
timeCAA = data{6};
timeCSA = data{7};
timeCSAoverCAA = data{8};

for kk = 1:numel(bytesCAA)
    if bytesCAA(kk) == 0 || bytesCSA(kk) == 0
        timeCSAoverCAA(kk) = NaN;
    end
end

if doStream; strTitle = 'Download times from Cluster Active/Science Archive, STREAMING';
else strTitle = 'Download times from Cluster Active/Science Archive'; end
plot(bytesCAA*1e-6,timeCSAoverCAA,'b*')
title(strTitle)
ylabel('t_{CSA}/t_{CAA}')
xlabel('file size [Mbytes]')
if doPrint; print('-dpng',[fileDir,startTimeStr,'_',streamStr]), end
end