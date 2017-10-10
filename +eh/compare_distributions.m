% Load data and make a variety of plots for electron distributions.
if 0
cd /Users/Cecilia/Data/BM/20070831
sclist=3:4;
%c_eval('peaDEFlux? = c_caa_distribution_data(''C?_CP_PEA_3DXPH_DEFlux'');',sclist)
%c_eval('peaDPFlux? = c_caa_distribution_data(''C?_CP_PEA_3DXPH_DPFlux'');',sclist)
%c_eval('peaPSD? = c_caa_distribution_data(''C?_CP_PEA_3DXPH_PSD'');',sclist)
%save('peaDEFlux.mat','peaDEFlux3','peaDEFlux4');
%save('peaDPFlux.mat','peaDPFlux3','peaDPFlux4');
%save('peaPSD.mat','peaPSD3','peaPSD4');
load peaDEFlux
load peaDPFlux
load peaPSD
load matlabE
end
%%

fig = figure(66); 
set(fig,'position',[560   209   562   739])
positionTop = [0.1300    0.4694    0.7750    0.4215];
positionBottom1 = [0.1300    0.2886    0.7750    0.1077];
positionBottom2 = [0.1300    0.1100    0.7750    0.1068];

h(1)=axes('position',positionTop);
h(2)=axes('position',positionBottom1);
h(3)=axes('position',positionBottom2);


tstart = toepoch([2007 08 31 10 17 35]);
tstop = toepoch([2007 08 31 10 17 50]); %tstart+20;peaPSD3.t(end);
dt = 2;

t1 = tstart; 
tint=cell(1,3);
numberOfTimeIntervals=0;
while t1+4*dt  < tstop
    numberOfTimeIntervals = numberOfTimeIntervals + 1;
    t2 = t1 + dt; tint{1}=[t1 t2];
    t3 = t2 + dt; tint{2}=[t2 t3];
    t4 = t3 + dt; tint{3}=[t3 t4];
    t5 = t4 + dt; tint{4}=[t4 t5];
    
    nint=4;
    tind=cell(nint,1);
    for k=1:nint
        [t,tind{k}] = irf_tlim(peaPSD3.t,tint{k});
        
        yy{k} = nanmean(squeeze(peaPSD3.p(tind{k},end,:)))';
        xx{k} = peaPSD3.en_cs;        
    end
    if 0
    loglog(h(1),xx{1},yy{1},'color',cn.colors(1,1),...
                xx{2},yy{2},'color',cn.colors(2,1),...
                xx{3},yy{3},'color',cn.colors(3,1),...
                xx{4},yy{4},'color',cn.colors(4,1))
    end
    loglog(h(1),xx{1},yy{1},...
                xx{2},yy{2},...
                xx{3},yy{3},...
                xx{4},yy{4})    
    
    title(h(1),num2str(numberOfTimeIntervals))
    t1=t3;
    
irf_plot(h(2),irf_tlim(diE3,[t1 t5]))
    irf_plot(h(3),irf_tlim(diE4,[t1 t5]))
    irf_zoom(h(2:3),'x',[t1 t5])
    irf_zoom(h(2:3),'y')
    c_eval('cn.print(''time_develeopment_?'')',numberOfTimeIntervals)
end
