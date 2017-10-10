function [ t0, corr ,corr_max ] = cn_xcorr(data1,data2,window,flag)
fs=1/(data1(2,1)-data1(1,1));
%fs=450;
ndata=size(data1,1);
switch flag
    case 'x'
        %size(data1(:,2))
        %size(data2(:,2))
        corr=xcorr(data1(:,2),data2(:,2),'coeff');
    case 'y'
        corr=xcorr(data1(:,3),data2(:,3),'coeff');
    case 'z'
        corr=xcorr(data1(:,4),data2(:,4),'coeff');
    case 'xy'
        corr=xcorr2(data1(:,2:3),data2(:,2:3));
        corr=corr(:,2);
end
    
OFF = fix(size(corr,1)/2)-window; % Have to take out a large number of values because corr looks like an M
ii = find(corr == max(corr(1+OFF:end-OFF)));

OFF=1; % Only use 3 points for polyfit
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
p = polyfit((ii-OFF:ii+OFF)',corr(ii-OFF:ii+OFF),2);
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
t0 = -p(2)/p(1)/2;

if abs(t0-ii)<1
    corr_max = polyval(p,t0);
else
    disp('Discarding parabolic fit!')
    t0 = ii;
    corr_max = corr(ii(1));
end

t0=t0-ndata;
t0=t0/fs;