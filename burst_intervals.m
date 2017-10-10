function burst_intervals(time,varargin)
% Calculate and plot distances between C3 and C4 from beginning of June to 
% end of October 

    sc_list=4; % Need only check one s/c since burst modes coincide
               % s/c 4 has no internal burst so no confusion
    coord_label='GSE'; % ?
    %%
        
switch size(varargin,1)
    case 0 % If only one time is given we find the next burst mode interval    
        indicator=1;
        
        while indicator 
        [tr,r] = caa_is_get('db.irfu.se:0', toepoch(time), 2, ic, 'ephemeris', 'position');
        c_eval('R?=[double(tr) double(r)''];',ic);
        clear tr r;
        
        
        %% Calculating the distance
        rx=R4(1,2)-R3(1,2);
        ry=R4(1,3)-R3(1,3);
        rz=R4(1,4)-R3(1,4);
        D34=sqrt(rx^2+ry^2+rz^2);
        disp(['D34 = ',num2str(D34),' km'])
        
        
    case 1  % If also a stop time is given we find all intervals wihtin this period
        stop_time=varargin{1,1};
        %%
        %time=[2007 06 10 10 10 10];
        %stop_time=[2007 11 10 10 10 10];
        
        %bursts=[0 0];
        
        burstomode0=0;
        start_time=time; 
        t=toepoch(time);
        tstop=toepoch(stop_time);
        
        t_vec=t:1:tstop;

        %% Looping over all times
        t0=t;
        t1=t0+100;
        dt=t1-t0;
        
        for n=1:length(t_vec)
            for ic=sc_list
                 [tr,r] = caa_is_get('db.irfu.se:0', t_vec(n), 1, ic, 'ephemeris', 'position');
                 %c_eval('R?=[double(tr) double(r)''];',ic);
                 %clear tr r;
            end
            t1=tr;
            f=1/(t1-t0);
            
            if f > 440
                burstmode1=1;
            else
                burstmode1=0;
            end
            
            if burstmode1>burstmode0 % we enter burst mode
                bursts
                
            elseif burstmode1<burstmode0 % we exit from burst mode
                
            else
            end
%%

            disp(['D34 = ',num2str(D34(n)),' km  Time = ', num2str(fromepoch(t_vec(n)))])
        end
        %% 
        D=[t_vec' D34'];
        %%
        h=irf_plot(2);
        irf_plot(h(2),D(:,[1 2]),'-*');
        ylabel('D_{34}')
end
