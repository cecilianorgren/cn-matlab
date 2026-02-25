function D34=dist34(time,stop_time)
% Calculate and plot distances between C3 and C4 from beginning of June to 
% end of October 

%time=[2007 06 10 10 10 10];
%stop_time=[2007 09 10 10 10 10];

sc_list=3:4;
coord_label='GSE';
start_time=time;
t=toepoch(time);
tstop=toepoch(stop_time);
%%
nt=120; % number of times
timestep=fix((tstop-t)/nt);
t_vec=t:timestep:tstop;

%% Looping over all times
D34=t_vec*0;

for n=1:length(t_vec) 
    
    for ic=sc_list
         [tr,r] = caa_is_get('db.irfu.se:0', t_vec(n), 1, ic, 'ephemeris', 'position');
         c_eval('R?=[double(tr) double(r)''];',ic);clear tr r;
    end
    
    rx=R4(1,2)-R3(1,2);
    ry=R4(1,3)-R3(1,3);
    rz=R4(1,4)-R3(1,4);
    
    D34(n)=sqrt(rx^2+ry^2+rz^2);
end

D=[t_vec' D34'];

%%
h=irf_plot(1);
irf_plot(h(1),D);
%plot(D34)
    
    
%disp(['D34 = ',num2str(D34),' km'])


