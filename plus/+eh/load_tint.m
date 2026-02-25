function [tint nscobs comments color] = load_tint
% EH.LOAD_TINT loads /Users/Cecilia/Research/EH/EH.txt' that contains time
%   intervals to electron phase space holes during the burst mode event on 
%   2008-08-31.
%
%   [tints, quality, comments]=EH.LOAD_TINT;
%       tints - cell array of time intervals in epoch
%       quality - 'quality level'
%           1 - Best
%           2 
%           3 
%           4 
%           5 
%           6 - Rubbish
%       comments - comments
%       color - 1x8 cell array of colors that can be used for marking
%               green, green, green, yellow, white, blue, red

[year1 month1 day1 hour1 min1 sec1 year2 month2 day2 hour2 min2 sec2 nscobs comments]...
    =textread('/Users/Cecilia/Research/EH/EH.txt','%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n]');
t11=[year1 month1 day1 hour1 min1 sec1];
t22=[year2 month2 day2 hour2 min2 sec2];
for k=1:size(year1,1); 
    tint{k,1}=[toepoch(t11(k,:)) toepoch(t22(k,:))];
end
color={'green','green','green','yellow','yellow','white','blue','red'};

