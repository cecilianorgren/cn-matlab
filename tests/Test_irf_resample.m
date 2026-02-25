%% Artificial data
utcT1 = '2002-03-04T09:30:00Z'; 
time = EpochTT(utcT1):1:(EpochTT(utcT1)+20);
f = 1.5/(time.stop-time.start);
data = repmat(cos(2*pi*f*(time-time.start)),1,3,3); %data = permute(data,[3 1 2]);


T = irf.ts_tensor_xyz(time,data);
V = irf.ts_vec_xyz(time,squeeze(data(:,:,1)));
S = irf.ts_scalar(time,squeeze(data(:,1,1)));

matV = [time.epochUnix squeeze(data(:,:,1))];

%%

h = irf_plot(5); isub = 1;

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,V)
irf_plot(hca,V.resample(((time):1:(time+2)),'max'),'ko')
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,V)
irf_plot(hca,V.resample(time(1:2:end)+0.5,'max'),'ko')
irf_plot(hca,V.resample(time(1:2:end)+0.5,'average'),'bo')
irf_plot(hca,V.resample(time(17)+0.1,'average'),'rs')
irf_plot(hca,V.resample(time.start:0.3:time.stop,'thresh'),'rs')
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,V)
irf_plot(hca,V.resample(time.start:0.3:time.stop,'thresh',5),'rs')
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,matV)
resampleTime = time.start:0.3:time.stop;
irf_plot(hca,irf_resamp(matV,resampleTime.epochUnix),'rs')
hold(hca,'off');
%%
intV = irf.ts_vec_xyz(time,int64(data(:,:,1)));
intV = intV.resample(time(1:2:end)); intV, class(intV.data)
intV = intV.resample(time.start:0.3:time.stop); intV, class(intV.data)
%%
intT = irf.ts_tensor_xyz(time,squeeze(int64(data)));
newT = intT.resample(time.start:0.3:time.stop), class(newT.data)
%% time,squeeze(data(:,:,1))

%% Real data
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
V = gseB1.tlim(tint);
T = gsePe1.tlim(tint);

time = V.time(fix(V.length/2));
h = irf_plot(5); isub = 1;

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,V)
irf_plot(hca,V.resample(((time):1:(time+3))+1),'ko')
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,V)
irf_plot(hca,V.resample(V.time(1:100:end)),'k.')
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,T.yy)
irf_plot(hca,T.yy.resample(T.time(1:100:end),'max'),'o')
irf_plot(hca,T.yy.resample(T.time(1:100:end),'average'),'o')
irf_plot(hca,T.yy.resample(T.time(1:100:end),'median'),'o')
irf_plot(hca,T.clone(T(1).time,-0.01).resample(T.time))
irf_legend(hca,{' ','max','average','median'},[0.98 0.95])
hold(hca,'off');

hca = h(isub); isub = isub + 1; hold(hca,'on');
irf_plot(hca,T.yy)
irf_plot(hca,T.resample(T.time(1:100:end),'max').yy,'o')
irf_plot(hca,T.yy.resample(T.time(1:100:end),'average'),'o')
irf_plot(hca,T.yy.resample(T.time(1:100:end),'median'),'o')
irf_plot(hca,T.clone(T(1).time,-0.01).resample(T.time))
irf_legend(hca,{' ','max','average','median'},[0.98 0.95])
hold(hca,'off');

irf_zoom(h,'x',tint)