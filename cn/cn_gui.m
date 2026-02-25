function cn_gui(t1,t2,E3_raw,E4_raw)

c_eval('E?=cn_toepoch(t1,t2,E?_raw);',3:4);

figure
h=irf_plot(2);
isub=1;

index=10;

if 1 % Allt
        h1=h(isub);isub=isub+1;
        irf_plot(h1,E3(:,[1 2]),'g'); hold(h1,'on');
        irf_plot(h1,E4(:,[1 2]),'b'); hold(h1,'on');
        irf_plot(h1,E3(index,[1 2]),'ko'); hold(h1,'on');
        ylabel(h1,'E_{X}[mV/m]');
        set(h1a,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(h1,{'C3','C4'},[0.02 0.05]); 
        set(h1,'ColorOrder',[0 0 0]);  
        irf_zoom(h1,'y');       
        ylim=get(h1,'ylim');
        xlim=get(h1,'xlim');
               
        h2=h(isub);isub=isub+1;
        irf_plot(h2,E3(:,[1 3]),'g'); hold(h2,'on');
        irf_plot(h2,E4(:,[1 3]),'b'); hold(h2,'on');
        irf_plot(h2,E3(index,[1 3]),'ko'); hold(h2,'on');
        
        ylabel(h2,'E_{Y,}[mV/m]');
        set(h2,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(h2,{'C3','C4'},[0.02 0.05]); 
        set(h2,'ColorOrder',[0 0 0]);   
        irf_zoom(h2,'y'); 
                
end
axpos=get(h2,'position')

uicontrol('Style', 'slider',...
        'Min',E3(1,1),'Max',E3(end,1),'Value',E3(1,1),...
        'Position', [65 20 495 20],...
        'Callback', {@time,[h1 h2]});
end

function time(hObj,event,ax) %#ok<INUSL>
    % Called to change time
    % when user moves the slider control 
    val = 51 - get(hObj,'Value');
    zlim(ax,[-val val]);
end