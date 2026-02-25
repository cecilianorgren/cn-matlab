% download msp list
urlList = 'http://jsoc1.bnsc.rl.ac.uk/msp/full_msp_ascii.lst';
msp = urlread(urlList);

% divide into cell arrays with burst mode intervals in text format
BM1s = regexp(msp,'(\d*)-(\d*)-(\d*)T(\d*):(\d*):(\d*)Z (\d*)-(\d*)-(\d*)T(\d*):(\d*):(\d*)Z B1','match');
BM2s = regexp(msp,'(\d*)-(\d*)-(\d*)T(\d*):(\d*):(\d*)Z (\d*)-(\d*)-(\d*)T(\d*):(\d*):(\d*)Z B2','match');

if 0 % write to textfile
    % BM1
    filePath = '/Users/Cecilia/Data/';
    fileName = 'BM1.txt';
    fileID = fopen([filePath fileName],'w');
    formatSpec= '%s\n';
    nBM1 = numel(BM1s);
    for k = 1:nBM1
        fprintf(fileID,formatSpec,BM1s{k}(1:41));
    end
    fclose(fileID);
    disp(' ')
    disp('List of BM1')
    disp('-----------------------------------------')
    type '/Users/Cecilia/Data/BM1.txt'

    % BM2
    filePath = '/Users/Cecilia/Data/';
    fileName = 'BM2.txt';
    fileID = fopen([filePath fileName],'w');
    formatSpec= '%s\n';
    nBM2 = numel(BM2s);
    for k = 1:nBM2
        fprintf(fileID,formatSpec,BM2s{k}(1:41));
    end
    fclose(fileID);
    disp(' ')
    disp('List of BM2')
    disp('-----------------------------------------')
    type '/Users/Cecilia/Data/BM2.txt'
end