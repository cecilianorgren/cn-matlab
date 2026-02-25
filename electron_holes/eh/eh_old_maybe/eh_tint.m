function [tint nscobs comments] = eh_tint
% EH_TINT loads /Users/Cecilia/Research/EH/EH.txt' that contains time
%   intervals to electron phase space holes during the burst mode event on 
%   2008-08-31.
%
%   [tints, quality, comments]=EH_TINT;
%       tints - cell array of time intervals in epoch
%       quality - 'quality level'
%           1 
%           2 
%           3 
%           4 
%           5 
%           6 

[year1 month1 day1 hour1 min1 sec1 year2 month2 day2 hour2 min2 sec2 nscobs comments]...
    =textread('/Users/Cecilia/Research/EH/EH.txt','%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n]');
t11=[year1 month1 day1 hour1 min1 sec1];
t22=[year2 month2 day2 hour2 min2 sec2];
for k=1:size(year1,1); 
    tint{k,1}=[toepoch(t11(k,:)) toepoch(t22(k,:))];
end

