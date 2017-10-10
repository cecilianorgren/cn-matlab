time=ud.tlim_mva;
time1=fromepoch(ud.tlim_mva(1));
time2=fromepoch(ud.tlim_mva(2));
vector=ud.v3;
l2l3=ud.l(2)/ud.l(3);
bminvar=ud.Xminvar;
bmin=sum(bminvar(:,4))/size(bminvar,1);
%%
myformat = '\n%4d %2d %2d %2d %2d %5.2f    %4d %2d %2d %2d %2d %5.2f    [%6.3f %6.3f %6.3f]    %5.1f    %5.1f';

% open the file with permission to append
fid = fopen('/Users/Cecilia/cn-matlab/nhat.txt','a');

% write values at end of file
fprintf(fid, myformat, [time1 time2 vector bmin l2l3]);

% close the file 
fclose(fid);
type nhat.txt

%%
time=ud.tlim_mva;
time1=fromepoch(time(1));
time2=fromepoch(time(2));
vector=ud.v3;
l=ud.l;


