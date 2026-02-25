function str = cn_time(epoch)
% transforms epoch to hhmmssmmm
%
%
time=fromepoch(epoch);

str=datestr(time,'HHMMSSFFF');