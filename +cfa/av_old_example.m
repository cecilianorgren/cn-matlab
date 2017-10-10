tmpDir = tempname;
mkdir(tmpDir);
cd(tmpDir);

if 1, % Example 1
    disp('_____ Example 1 ______________________')
    disp('==== CAA ====')
    % ********** CAA **************
    %tic;caa_download('2001-02-01T00:00:00.000Z/2001-02-01T02:00:00.000Z','C3_CP_RAP_I3D_H','nowildcard'); toc
    % last value: Elapsed time is 7.855032 seconds.
    % ********** CFA **************
    disp('==== CFA ====')
    %temp_file=tempname;
    temp_file='delme';
    url_line='http://cfaint.esac.esa.int/cfa/aio/product-action?DATASET_ID=C3_CP_RAP_I3D_H&START_DATE=2001-02-01T00:00:00Z&END_DATE=2001-02-01T02:00:00Z';
    disp(url_line);
    tic; urlwrite(url_line,temp_file);
    %delete(temp_file);
    toc
    % last value: Elapsed time is 27.422193 seconds.
end
ls -l
%rmdir(tmpDir,'s');