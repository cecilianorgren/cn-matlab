function D = dist34(time,varargin)
% Calculate and plot distances between C3 and C4 from beginning of June to 
% end of October 

    sc_list=3:4;
    coord_label='GSE';
    %%
    
switch size(varargin,1)
    case 0 % If only one time is given
        burstmode=0;
        for ic=sc_list % Looping over C3 and C4
             [tr,r] = caa_is_get('db.irfu.se:0', toepoch(time), 2, ic, 'ephemeris', 'position');
             c_eval('R?=[double(tr) double(r)''];',ic);
             clear tr r;
        end
        %% Calculating the distance
        rx=R4(1,2)-R3(1,2);
        ry=R4(1,3)-R3(1,3);
        rz=R4(1,4)-R3(1,4);
        D34=sqrt(rx^2+ry^2+rz^2);
        disp(['D34 = ',num2str(D34),' km'])
        D=D34;
    case 1  % If also a stop time is given   
        stop_time=varargin{1,1};
        %%
        %time=[2007 06 10 10 10 10];
        %stop_time=[2007 11 10 10 10 10];
        
        start_time=time; 
        t=toepoch(time);
        tstop=toepoch(stop_time);
        
        nt=420; % number of times
        timestep=fix((tstop-t)/nt);
        t_vec=t:timestep:tstop;

        %% Looping over all times
        D34=t_vec*0;

        for n=1:length(t_vec) 
            for ic=sc_list
                 [tr,r] = caa_is_get('db.irfu.se:0', t_vec(n), 1, ic, 'ephemeris', 'position');
                 c_eval('R?=[double(tr) double(r)''];',ic);
                 clear tr r;
            end 
%%
            rx=R4(1,2)-R3(1,2);
            ry=R4(1,3)-R3(1,3);
            rz=R4(1,4)-R3(1,4);
            D34(n)=sqrt(rx^2+ry^2+rz^2);
            disp(['D34 = ',num2str(D34(n)),' km  Time = ', num2str(fromepoch(t_vec(n)))])
        end
        %% 
        D=[t_vec' D34'];
        %%
          h=irf_plot(2);
          irf_plot(h(2),D(:,[1 2]),'-*');
          ylabel('D_{34} [km]')
end