%% With time series
c_eval('R? = mvaR?.resample(mvaB1);',1:4)

c_eval('B? = mvaB?.resample(mvaB1);',1:4)
[curvB,BB]=c_4_grad('R?','B?','curvature'); curvB.name = 'curv B';
curvBradius = 1/curvB.abs; curvBradius.name = 'R_c';

c_eval('B?y0 = B?; B?y0.data(:,2) = B?.data(:,2)*0;',1:4)
[curvBy0,BBy0]=c_4_grad('R?','B?y0','curvature'); curvBy0.name = 'curv B (By=0)';
curvBradiusy0 = 1/curvBy0.abs; curvBradiusy0.name = 'R_c';
c_eval('B?noguide = B?; B?noguide.data(:,2) = B?.data(:,2)-4;',1:4)
[curvBnoguide,BBnoguide]=c_4_grad('R?','B?noguide','curvature'); curvBnoguide.name = 'curv B (noguide)';
curvBradiusnoguide = 1/curvBnoguide.abs; curvBradiusnoguide.name = 'R_c';


h = irf_plot({B1,B1.abs,curvB,curvB.abs,curvBy0,curvBy0.abs,B1noguide,B1noguide.abs,curvBnoguide,curvBnoguide.abs});

%%
h = irf_plot({B1,B2,B3,B4,curvB,curvB.abs,mvaJcurl});
%h(4).YLim = [0 100];
%h(5).YLim = [0 10];
%[gradBdrift,BB]=c_4_grad('R?','B?','drift_grad');
%[curvBdrift,BB]=c_4_grad('R?','B?','drift_curl');

%h = irf_plot({mvaB1,curvB,curvB.abs,1/curvB.abs,1/curvB.abs/re1,gradBdrift,curvBdrift});