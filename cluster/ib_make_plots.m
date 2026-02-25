% ib_make_plots

% Internal burst data info:
% column 1: time
% column 2: sc pot
% column 2: sc pot
% column 2: sc pot

cd /Users/Cecilia/Data/IB;
load date_ib;
load tint_ib;
sclist=3:4;

for k=21;30:length(str3)
    eval(['cd /Users/Cecilia/Data/IB/',str3{k}]);
    cd
    
    % load data
    c_eval('EB?=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sclist);
    c_eval('B?=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);
    
    if 1
    % plot figure
    h=irf_plot(4); isub=1;
    set(gcf,'PaperPositionMode','auto');
    title(h(1),str3{k})
    if 0 % Plot all, E3, E4 in 2 panels total.
    if ~isempty(EB3)
        hca=h(isub);isub=isub+1;
        irf_plot(hca,EB3(:,[1 2:end]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if ~isempty(EB4)
        hca=h(isub);isub=isub+1;
        irf_plot(hca,EB4(:,[1 2:end]));
        ylabel(hca,'C4');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    irf_zoom()
    end    
    if 0 % plot all separately, 8 panels total
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 2]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 3]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 4]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 5]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 6]));
        ylabel(hca,'C3')
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 7]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 8]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3(:,[1 9]));
        ylabel(hca,'C3');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    end   
    if 0
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E4(:,[1 2]));
        ylabel(hca,'C4');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 0
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E4(:,[1 3]));
        ylabel(hca,'C4');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 0
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E4(:,[1 4]));
        ylabel(hca,'C4');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1 % Compare E between C3 and C4
        if 1
            hca=h(isub);isub=isub+1;
            irf_plot(hca,B3(:,[1:4]));
            ylabel(hca,'B C3 [nT]');
            irf_legend(hca,{'x','y','z'},[0.95 0.95]);
        end
        if 1
            hca=h(isub);isub=isub+1;
            irf_plot(hca,B4(:,[1:4]));
            ylabel(hca,'B C4 [nT]');
            irf_legend(hca,{'x','y','z'},[0.95 0.95]);
        end
        if 1
            hca=h(isub);isub=isub+1;
            c_pl_tx(hca,EB3,EB4,EB3,EB4,2);
            ylabel(hca,'E_x [mV/m]');
            irf_legend(hca,{'C3','C4'},[0.95 0.95]);
        end
        if 1
            hca=h(isub);isub=isub+1;
            c_pl_tx(hca,EB3,EB4,EB3,EB4,3);
            ylabel(hca,'E_y [mV/m]');
            irf_legend(hca,{'C3','C4'},[0.95 0.95]);
        end
    irf_zoom(h,'x',[min([EB3(:,1);EB4(:,1)])-2 max([EB3(:,1);EB4(:,1)])+2])    
    end
    
    print -dpng E&B.png
    %close gcf
    end
end