function savet(str)
figdata=get(gcf,'userdata');
interv=get(gca,'xlim');
time1=figdata.t_start_epoch+interv(1);time1=fromepoch(time1);
time2=figdata.t_start_epoch+interv(2);time2=fromepoch(time2);

% open the file with permission to append
fid = fopen('/Users/Cecilia/cn-matlab/tlim.txt','a');

% write values at end of file
myformat = 't1=[%4d %2d %2d %2d %2d %7.4f];   t2=[%4d %2d %2d %2d %2d %7.4f];';
fprintf(fid, myformat, [time1 time2]);
fprintf(fid, '  ''%s''; \n', str);

% close the file 
fclose(fid);
type tlim.txt