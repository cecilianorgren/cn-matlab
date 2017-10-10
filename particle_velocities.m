h2=irf_plot(14);
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
            ylabel(hca,'B to-\newline spinplane \newline angle'); grid(hca,'off');
            irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,codVi3);hold(hca,'on');
            ylabel(hca,'Codif3 \newline V_i GSE'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,codVi4);hold(hca,'on');
            ylabel(hca,'Codif4 \newline V_i GSE'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,hiaVi3);hold(hca,'on');
            ylabel(hca,'Hia3 \newline V_i GSE'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,gsmExB3);hold(hca,'on');
            ylabel(hca,'V_{ExB} \newline GSM C3'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,gsmExB4);hold(hca,'on');
            ylabel(hca,'V_{ExB} \newline GSM C4'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,peaVe3);hold(hca,'on');
            ylabel(hca,'PEA V_{e} \newline GSM? C3'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 1 % 
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,peaVe4);hold(hca,'on');
            ylabel(hca,'PEA V_{e} \newline GSM? C4'); grid(hca,'off');
            irf_legend(hca,{'x','y','z'},[0.02 0.05]);
        end
        if 0
            hca=h2(isub); isub=isub+1;
            irf_plot(hca,spec{4},'colorbarlabel',spec{4}.p_units);
            set(hca,'yscale','log','ytick',[1e0 1e1 1e2 1e3 1e4])
            ylabel(hca,['Energy [',spec{4}.f_units,']'])
            caxis(hca,[-2 2])
            hc=colorbar('peer',hca);
            ylabel(hc,spec{4}.p_units)
            colormap(hca,xcm)
        end
        if 0
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
        irf_pl_mark(h,tint,'red')
        linkaxes([h2 hbeta],'x');
    end