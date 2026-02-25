function tint = load_tint_ov
% EH.LOAD_TINT_OV Returns time interval for event overview in epoch.
%   2007-08-31 10:14:30 - 2007-08-31 10:14:30
%
%   tint = EH.LOAD_TINT_OV;

t1=[2007 08 31 10 14 30];
t2=[2007 08 31 10 19 30];
tint=toepoch([t1;t2])';