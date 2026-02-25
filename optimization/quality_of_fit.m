function [Q_xcorr,Q_area] = quality_of_fit(A1,A2)
% Gives quality factor Q of data series A1 and A2.
% Q is a combination of waveform correlation and the squared area between series.

Q_xcorr = xcorr(A1,A2,'coeff',0);
QA_area_std = std(A1);
Q_area = sum(sqrt((A1-mean(A1)-(A2-mean(A2))).^2))/numel(A1);
%Q_area = QA_area/QA_area_std;