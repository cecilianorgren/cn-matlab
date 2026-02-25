%% Figure with zig-zag in 4s intervals + overview
% load data
%c_eval('res?=c_caa_construct_subspin_res_data(irf_ssub(''Data__C?_CP_PEA_3DXPH_PSD'',?));',3:4);
load matlabE;
load matlabB;
load matlabBeta;
load matlabN;
varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',4);
[var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
DEF=squeeze(nansum(nansum(varmat.data,2),3));
energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);

c_eval('res?=c_caa_construct_subspin_res_data(irf_ssub(''Data__C?_CP_PEA_3DXPH_PSD'',?));',3:4);


c_eval('diB?=c_coord_trans(''gsm'',''dsi'',gsmB?fgm,''cl_id'',?);',3:4);
c_eval('[facE?,angle?]=irf_edb(diE?,diB?,0);',3:4);

%% utg? fr?n irf_plot
% custom colormap
it=0:.02:1;it=it';
xcm=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 ...
    flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]];
clear it;

clear h h1 h2 h3 hca p k
% define time step etc.
tint_ov=[toepoch([2007 08 31 10 13 30]) toepoch([2007 08 31 10 19 30])];
[tint_eh, scncobs, ~]=eh_tint;
dt=20; % s
t1=tint_ov(1);
t2=t1+dt;
mmm=1;
while t1<tint_ov(end)
    mmm=0;
% set up figure
figure('position',[1 -21 1168 955])
np=18;
wid=0.37; mid=0.55;
left=0.08;
h1=irf_plot(np);
h1pos=get(h1,'position');
for k=1:np
    h1pos{k}(3)=wid;
    h1pos{k}(1)=left;
    h2pos{k}=h1pos{k};
    h2pos{k}(1)=mid;
end
for k=1:np
    set(h1(k),'position',h1pos{k});
    if k<10
        h2(k)=axes('position',h2pos{k});
    elseif k>10   
        %k-10
        h3(k-10)=axes('position',h2pos{k});
    end
end



% make the plot
% c3 has less energy levels, c3.en(1)=c4.en(5)
enindex=[2 5 8 11 14 17 20 23 26]; % for c3 between 1 and 26
ne=length(enindex);
ic=3:4;

    tint=[t1 t2]; t1=t1+dt; t2=t2+dt;
    if 1,   % PANEL: PEACE 3DXPH_DEFlux high res angular spectrogra,
        for l=1:2;           
            c_eval('res=res?;',ic(l));
            [delmett,ind]=irf_tlim(res.tt,tint);
            specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);        
            for p=1:ne; % pitch angle spectrogram for given energy
                if l==1 % c3
                    eind=enindex(end-(p-1));
                else % c4
                    eind=enindex(end-(p-1))+4;
                end
                eind;
                hca=h1((p-1)*2+l);
                specrec.f=res.theta;specrec.f_label='PA';
                specrec.p=res.pitch_angle(ind,:);
                specrec.f_label=[specrec.f_label '  \newline' num2str(res.en(eind),4) 'eV'];
                specrec.p=log10(res.data(ind,:,eind));
                hs((p-1)*2+l)=irf_spectrogram(hca,specrec);
                set(hca,'ytick',[30 60 90 120 150]);
                %caxis(hca,[-2.5 0.5]);
                irf_legend(hca,['C' num2str(ic(l))],[0.98,0.98]);
            end
        end
        for p=2:2:2*ne
            cmins=min([caxis(hs(p-1));caxis(hs(p))]);
            cmaxs=max([caxis(hs(p-1));caxis(hs(p))]);
            clims=[min(cmins) max(cmaxs)];
            caxis(hs(p),clims);
            caxis(hs(p-1),clims);
        end
        
        title(h1(1),'Data__C?_CP_PEA_3DXPH_PSD')
        irf_timeaxis(h1(end))
        linkaxes(h1,'x')
    %isub=max(size(enindex))+1;
    end
    if 1 % overview plot on the RHS
        isub=1;
        % overview with markings
        if 1 % gsmB4
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,gsmB4fgm);
            ylabel(hca,'B [nT] \newline GSM C4');
            irf_legend(gca,{'x','y','z'},[0.98 0.95]);
            set(hca,'ylim',[-5 35]); grid(hca,'off')
        end
        if 1 % Electron DEF
            hca=h2(isub); isub=isub+1;
            ud = get(gcf,'userdata');
            t_st_e = double(ud.t_start_epoch);
            pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
            set(hca,'yscale','log','tickdir','out');
            set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
            ylabel(hca,'E_e [eV]');
            hcb = colorbar('peer',hca);
            irf_colormap(hcb,'space')
            ylabel(hcb,{'log_{10} dEF';varunits},'fontsize',14)
            shading(hca,'flat')
            caxis(hca,[6.5 10.3]);
            hcaylim=get(hca,'ylim');
            set(hca,'ylim',[60 hcaylim(2)])
            irf_legend(hca,'C4',[0.98 0.05],'color','k')  
        end
        if 1 %x Peace electron density and plasma beta (1 panel)    
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,peaNe4);
            set(hca,'box','off')
            hca(2) = axes('Position',get(hca(1),'Position'));
            hbeta=hca(2);
            irf_plot(hca(2),betaEH,'r')
            set(hca(2),'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
            set(hca(2),'YAxisLocation','right','YScale','log',...
                'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
            ylabel(hca(2),'\beta');
            set(hca(2),'Color','none','box','off'); % color of axis
            set(hca(2),'XColor','k','YColor','k'); % color of axis lines and numbers

            %h2=hca(2);

            irf_timeaxis(hca(2),'nolabels')
            irf_timeaxis(hca(1),'nolabels')

            irf_legend(hca(1),'n_e',[0.62 0.15],'color','k')
            irf_legend(hca(2),'\beta',[0.5 0.62],'color','r')
            ylabel(hca(1),'n_e [cm^{-3}]');
            grid(hca(1),'off');grid(hca(2),'off')
            set(hca(2),'ylim',[0.001 max(betaEH(:,2))*1.1])
            linkaxes(hca,'x')
            %set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
            %irf_legend(hca,{'PEACE'},[0.02 0.2]);

        end
        if 1 % E C3 ISR2 (1 panel)
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,irf_tlim(diE3(:,[1 2 3]),tint_ov(1)-1,tint_ov(2)+1));hold(hca,'on');
            ylabel(hca,'E [mV/m] \newline ISR2 C3');grid(hca,'off');
            irf_legend(hca,{'x','y'},[0.02 0.05]);
        end
        if 1 % E C4 ISR2 (1 panel)
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,irf_tlim(diE4(:,[1 2 3]),tint_ov(1)-1,tint_ov(2)+1));hold(hca,'on');
            ylabel(hca,'E [mV/m] \newline ISR2 C4');grid(hca,'off');
            irf_legend(hca,{'x','y'},[0.02 0.05]);
        end

        if 1 % E C4C4 Fac_z (1 panel)
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,{irf_tlim(facE3(:,[1 4]),tint_ov(1)-1,tint_ov(2)+1),...
                          irf_tlim(facE4(:,[1 4]),tint_ov(1)-1,tint_ov(2)+1)},'comp');hold(hca,'on');
            ylabel(hca,'E [mV/m]\newline FAC Z');grid(hca,'off');
            irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        end
        if 1 % E C4C4 BtoSpinplane angle (1 panel)
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,{irf_tlim([facE3(:,1) angle3],tint_ov(1)-1,tint_ov(2)+1),...
                          irf_tlim([facE4(:,1) angle4],tint_ov(1)-1,tint_ov(2)+1)},'comp');hold(hca,'on');
            ylabel(hca,'B to spin-\newline plane angle'); grid(hca,'off');
            irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        end
        if 1
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,spec{4},'colorbarlabel',spec{4}.p_units);
            set(hca,'yscale','log','ytick',[1e0 1e1 1e2 1e3 1e4])
            ylabel(hca,['Energy [',spec{4}.f_units,']'])
            caxis(hca,[-2 2])
            hc=colorbar('peer',hca);
            ylabel(hc,spec{4}.p_units)
            colormap(hca,xcm)
        end
        if 1
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,spec{5},'colorbarlabel',spec{5}.p_units);
            set(hca,'yscale','log','ytick',[1e0 1e1 1e2 1e3 1e4])
            ylabel(hca,['Energy [',spec{5}.f_units,']'])
            caxis(hca,[-2 2])
            hc=colorbar('peer',hca);
            ylabel(hc,spec{5}.p_units)
            colormap(hca,xcm)
        end
        irf_timeaxis(h2)
        irf_zoom([h2 hbeta],'x',tint_ov)
        cn_plot_axis_align([h2 hbeta]);
        irf_pl_mark(h2,tint,'red')
        linkaxes([h2 hbeta],'x');
    end
    % add some zoomed in quantities
    if 1
        isub=1;
        if 1 %E
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,{irf_tlim(diE3(:,[1 2]),tint),irf_tlim(diE4(:,[1 2]),tint)},'comp');
            irf_legend(hca,{'C3','C4'},[0.02 0.9]); ylabel(hca,'E_{X,ISR2}');grid(hca,'off');
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,{irf_tlim(diE3(:,[1 3]),tint),irf_tlim(diE4(:,[1 3]),tint)},'comp');
            irf_legend(hca,{'C3','C4'},[0.02 0.9]); ylabel(hca,'E_{Y,ISR2} ');grid(hca,'off');
        end
        if 1 % gsmB
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,{irf_tlim(gsmB3(:,[1 2]),tint),irf_tlim(gsmB4(:,[1 2]),tint)},'comp');
            irf_legend(hca,{'C3','C4'},[0.02 0.9]); 
            ylabel(hca,'B_{X,GSM}');grid(hca,'off');
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,{irf_tlim(gsmB3(:,[1 3]),tint),irf_tlim(gsmB4(:,[1 3]),tint)},'comp');
            irf_legend(hca,{'C3','C4'},[0.02 0.9]); 
            ylabel(hca,'B_{Y,GSM}');grid(hca,'off');
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,{irf_tlim(gsmB3(:,[1 4]),tint),irf_tlim(gsmB4(:,[1 4]),tint)},'comp');
            irf_legend(hca,{'C3','C4'},[0.02 0.9]); 
            ylabel(hca,'B_{Z,GSM}');grid(hca,'off');        
        end
        if 1 % deltaB
            hca=h3(isub);isub=isub+1; 
            irf_plot(hca,irf_add(1,irf_tlim(gsmB4,tint),-1,irf_tlim(gsmB3,tint)));
            irf_legend(hca,{'x','y','z'},[0.02 0.9]); 
            ylabel(hca,'B4-B3_{GSM}');grid(hca,'off');              
        end
            irf_zoom([h3],'x',tint)
            irf_zoom(h3,'y')
            linkaxes(h3,'x')
    end
    if 1
        % mark all eh's, and especially the presented one.
        % prio 1/2/3/4/5/rubbish/dl/presented interval
        color={'green','green','green','yellow','yellow','white','blue','red'};
        for p = 1:size(tint_eh,1)        
            cin=scncobs(p);        
            irf_pl_mark([h],tint_eh{p},color{cin})               
        end    
    end
    set(gcf,'PaperPositionMode','auto');
    eval(['print -dpng /Users/Cecilia/EH/Pics/PAD_OV_20s_mark__anis_T',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(2)),'HHMMSSFFF'),'.png'])
    close gcf
end