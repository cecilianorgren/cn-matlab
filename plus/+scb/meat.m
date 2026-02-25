%slCharacterEncoding('UTF-8')
% ? - char(229)
% ? - char(228)
% ? - char(246)
% ? - char(197)
% ? - char(196)
% ? - char(214)

set(0,'defaultAxesFontSize',14);
set(0,'defaultTextFontSize',14);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');
%% Slaktade kycklingar
fid    = fopen('/Users/Cecilia/Statistik/Slaktadekycklingar.txt');
format = '%*s%*s%f%f';
data   = textscan(fid,format,'headerlines',1);
fclose(fid);
plot(data{1},data{2}/1000);
set(gca,'ylim',[0 10^2])
title('Antal slaktade kycklingar, miljoner ')
xlabel([char(197) 'r'])
ylabel('Antal kycklingar [miljoner] ')
%% Slaktade kycklingar
fid    = fopen('/Users/Cecilia/Statistik/tonfjaderkott.txt');
format = '%*s%f%f';
data   = textscan(fid,format,'headerlines',1);
fclose(fid);
plot(data{1},data{2}/1000);
%set(gca,'ylim',[0 10^2])
title(['Produktion av fj' char(228) 'derk' char(246) 'tt, 1000 ton'])
xlabel([char(197) 'r'])
ylabel('Produktion [1000 ton]')
%% Produktion av k?tt, 2 slag
fid    = fopen('/Users/Cecilia/Statistik/prodavkott2_1000ton.txt');
format = '%*s%*s%*s%*s%*s%*s%*s%*s%f-%f%f';
data   = textscan(fid,format,9,'headerlines',1);
year = round((data{2}+data{1})/2);
beefhorsesheeplamb = data{3}; 
%%
format = '%*s%f-%f%f';
data   = textscan(fid,format,9,'headerlines',1);
pork = data{3};
%%
fclose(fid);
plot(year,beefhorsesheeplamb,year,pork);
hl=legend(['N' char(246) 'tkreatur, h' char(228) 'st, f' char(229) 'r och lamm'],'Gris')
set(hl,'edgecolor',[1 1 1],'location','southeast','box','off')
%set(gca,'ylim',[0 10^2])
title(['Produktion av k' char(246) 'tt, 1000 ton'])
xlabel([char(197) 'r'])
ylabel('Produktion [1000 ton]')
%% Produktion av k?tt, 5 slag
fid    = fopen('/Users/Cecilia/Statistik/prodavkott1000ton.txt');
format = '%*s%f/%*f%f';
data   = textscan(fid,format,65,'headerlines',1);
ar = data{1};
storboskap = data{2};
format = '%*s%f/%*f%f';
data   = textscan(fid,format,65,'headerlines',1);
kalv = data{2};
format = '%*s%*s%*s%f/%*f%f';
data   = textscan(fid,format,65,'headerlines',1);
hast = data{2};
format = '%*s%*s%*s%f/%*f%f';
data   = textscan(fid,format,65,'headerlines',1);
far = data{2};
format = '%*s%f/%*f%f';
data   = textscan(fid,format,65,'headerlines',1);
svin = data{2};
fclose(fid);
plot(ar,storboskap,ar,kalv,ar,hast,ar,far,ar,svin);
hl=legend('Storboskap','Kalv',['H' char(228) 'st inkl. f' char(246) 'l'],['F' char(229) 'r inkl. lamm'],'Svin')
set(hl,'edgecolor',[1 1 1],'location','northwest','box','off')
title(['Produktion av k' char(246) 'tt, 1000 ton'])
xlabel([char(197) 'r'])
ylabel('Produktion [1000 ton]')
%% Kombination av 2 slag och 5 slag
plot(year,beefhorsesheeplamb,'b',year,pork,'g'); hold(gca,'on')
plot(ar,storboskap+kalv+hast+far,ar,svin,'g'); hold(gca,'off')
%legend('Storboskap','Kalv','Hast inkl. fol','Far inkl. lamm','Svin')
leg_a1 = {['N' char(246) 'tkreatur'],['H' char(228) 'st'],['F' char(229) 'r'],'Lamm'};
leg_a2 = {'Storboskap','Kalv',['H' char(228) 'st inkl. f' char(246) 'l'],['F' char(229) 'r inkl. lamm'],'Svin'};
leg_b1 = {'Gris'};
leg_b2 = {'Svin'};
text(1910,50,leg_a1)
text(1930,270,leg_a2)
text(1880,120,leg_b1)
text(1975,110,leg_b2)
title(['Produktion av k' char(246) 'tt, 1000 ton'])
xlabel([char(197) 'r'])
ylabel('Produktion [1000 ton]')


