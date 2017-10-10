%% corrcoef
corcoef.x=corrcoef(bE3resamp(1:end-1,2),bE4(2:end,2)); % x
corcoef.y=corrcoef(bE3resamp(1:end-1,3),bE4(2:end,3)); % y
corcoef.x_det=corrcoef(detrend(bE3resamp(1:end-1,2),'constant'),detrend(bE4(2:end,2),'constant')); % x
corcoef.y_det=corrcoef(detrend(bE3resamp(1:end-1,3),'constant'),detrend(bE4(2:end,3),'constant')); % y

%% corrcoef to xcorr
[M,P] = size(bE3resamp(1:end-1,2));
c=xcov(bE3resamp(1:end-1,2),bE4(2:end,2),'coef');
c0 = zeros(P); c0(:) = c(M,:)
%% xcorr
[xcorrr.x xcorrr.x_lags]=xcorr(bE3resamp(:,2),bE4(:,2),'coeff'); % x
[xcorrr.y xcorrr.y_lags]=xcorr(bE3resamp(:,3),bE4(:,3),'coeff'); % y
[xcorrr.x_det xcorrr.x_det_lags]=xcorr(detrend(bE3resamp(1:end-1,2),'constant'),detrend(bE4(2:end,2),'constant'),'coeff'); % x
[xcorrr.y_det xcorrr.y_det_lags]=xcorr(detrend(bE3resamp(1:end-1,3),'constant'),detrend(bE4(2:end,3),'constant'),'coeff'); % y

figure;plot(xcorrr.x_det_lags,xcorrr.x_det,xcorrr.y_det_lags,xcorrr.y_det);legend('x','y');title('Correlation coefficients, detrended');xlabel('Correlation coefficients');ylabel('lag')
figure;plot(xcorrr.x_lags,xcorrr.x,xcorrr.y_lags,xcorrr.y);legend('x','y');title('Correlation coefficients');xlabel('Correlation coefficients');ylabel('lag')

%% try own functions to see how they are correlated
x = 0:0.01:10;
X = sin(x);
Y = cos(x); % 90 degree phase
Y1= -sin(x); % 180 phase + higher amplitude -1
Y2= -sin(x)*2; % 180 phase + higher amplitude -1
Y3= -sin(x)+2; % 180 phase + DC component -1
corrcoef(X,Y)

%% xcorr matlab example
x=0:0.01:10;
X = sin(x);
Y=cos(x);
xcor=xcorr(X,Y,'coeff');  
max(r)

%% xcorr detrend
xcorr.xdet=xcorr(detrend(bE3resamp(:,2)),detrend(bE4(:,2)),'coeff'); % x
xcorr.ydet=xcorr(detrend(bE3resamp(:,3)),detrend(bE4(:,3)),'coeff'); % y

%% xcov
xcov.x=xcov(bE3resamp(:,2),bE4(:,2)); % x
xcov.y=xcov(bE3resamp(:,3),bE4(:,3)); % y
