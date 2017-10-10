fid    = fopen('/Users/Cecilia/Statistik/Slaktadekycklingar.txt');
format = '%*s%*s%f%f';
data   = textscan(fid,format,'headerlines',1)
%%
h = plot(data{1},data{2}/1000);
%ymax = ceil(max(data{2}));
set(gca,'ylim',[0 10^2])%,'xtick',[min(data{1}) max(data{1})]);