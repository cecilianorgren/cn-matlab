function plot_other(time1,time2,time,unit,gseB3,gseB4)
switch unit
    case 's'
        time1(6)=time1(6)-time;
        time2(6)=time2(6)+time;
        t1=time1;
        t2=time2;
    case 'm'        
        time1(5)=time1(5)-time;
        time2(5)=time2(5)+time;
        t1=time1;
        t2=time2;
end

B3=cn_toepoch(t1,t2,gseB3);
B4=cn_toepoch(t1,t2,gseB4);
figure; 
h=irf_plot(4);

%irf_plot(h(1),(B3))
%irf_plot(h(2),(B4))
irf_plot({B3,B4},'comp');
irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)])
irf_plot_axis_align
irf_timeaxis(h,'usefig');
title(h(1),'B [nT] (x,y,z) C3 blue C4 green')