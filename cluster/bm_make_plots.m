% bm_make_plots

[t1,t2,tint_comments]=textread('/Users/Cecilia/Exjobb/BM/BM.txt','%s%s%[^\n]');
for j=1:size(t1,1)
    tint{j}=[iso2epoch(t1{j}) iso2epoch(t1{j})];
    datestring{j}=datestr(fromepoch(tint{j}(1)),'yyyymmdd');
end
n_tint=j;
clear t1 t2 tint_comments j
%%
cd /Users/Cecilia/Data/BM;
sclist=3:4;
%%
for k=1;    
    eval(['cd /Users/Cecilia/Data/BM/',datestring{k}]);
    % ls caa
    
    % load data
    if ~exist('EB.mat','file')
        c_eval('gseB?=c_caa_var_get(''B_vec_xyz_gse__C3_CP_FGM_FULL'',''mat'');',sclist);
        c_eval('c_coord_trans(''GSE'',''ISR2'',gseB?,''cl_id'',?);',sclist)
        c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
        c_eval('diB?=irf_abs(diB?);',sclist);
        save('EB','diE3','diE4','diB3','diB4');
    elseif ~exist('diE3','var')
        load EB;
    end
    
    if 1 % plot overview figure with B, E
    h=irf_plot(5,'newfigure'); isub=1;
    setupfigure;
    
    % Compare E between C3 and C4
    if 1 % B3
        hca=h(isub);isub=isub+1;
        irf_plot(hca,diB3);
        ylabel(hca,'B C3 [nT]');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1 % B4
        hca=h(isub);isub=isub+1;
        irf_plot(hca,diB4);
        ylabel(hca,'B C4 [nT]');
        irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    end
    if 1 % Ex
        hca=h(isub);isub=isub+1;
        irf_plot(hca,{diE3(:,[1 2]),diE4(:,[1 2])},'comp');
        ylabel(hca,'E_x [mV/m]');
        irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    end
    if 1 % Ey
        hca=h(isub);isub=isub+1;
        irf_plot(hca,{diE3(:,[1 3]),diE4(:,[1 3])},'comp');
        ylabel(hca,'E_y [mV/m]');
        irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    end    
    if 1 % Ez
        hca=h(isub);isub=isub+1;
        irf_plot(hca,{diE3(:,[1 4]),diE4(:,[1 4])},'comp');
        ylabel(hca,'E_z [mV/m]');
        irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    end
    title(h(1),[datestring{k},' ISR2'])  
    irf_zoom(h,'x',[diE3(1,1) diE3(end,1)])    
    print -dpng E&B_ov.png
    end
    
    if 1 % print zoomed in plots
    dt=60*3;
    t1=diE3(1,1); t2=t1+dt;
    while t1<diE3(end,1);        
        irf_zoom(h,'x',[t1 t2]);
        irf_zoom(h,'y')
        eval(['print -dpng E&B_',datestr([fromepoch(t1)],'HHMMSS'),'.png'])
        t1=t2; t2=t2+dt;    
    end
    end
    close gcf
end