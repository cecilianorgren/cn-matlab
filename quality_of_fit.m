function quality_factor = quality_of_fit(A1,A2)
% Gives quality factor Q of data series A1 and A2.
% Q is a combination of waveform correlation and the squared area between series.

QA_xcorr = xcorr(A1,A2,'coeff',0);
QA_area_std = std(A1);
QA_area = sqrt((A1-A2).^2);
QA_area = QA_area/QA_area_std;

quality_factor = QA_xcorr;