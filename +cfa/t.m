function out = t(t)
% CFA.TIME Transform time form isdat epoch to iso format.
%   timeInIso = CFA.TIME(timeInEpoch);

tintiso = irf_time([t t],'tint2iso');
nChar=numel(tintiso);
out = tintiso(1:floor(nChar/2));

