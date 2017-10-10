%t1=[2007 08 31 10 17 45 50]; t2=[2007 08 31 10 17 45 90];
t1=[2007 08 31 10 17 46 20]; t2=[2007 08 31 10 17 46 70];
t1=[2007 08 31 10 17 46 46]; t2=[2007 08 31 10 17 46 60];
t1=[2007 08 31 10 17 45 50]; t2=[2007 08 31 10 17 46 70];
t1=[2007 08 31 10 17 45 50]; t2=[2007 08 31 10 17 45 90]; % first spike

%%
err=0.1;
c_eval('e?=cn_toepoch(t1,t2,diE?);',3:4)
c_eval('b?=cn_toepoch(t1,t2,gseB?);',3:4)
c_eval('diB?=c_coord_trans(''gse'',''dsi'',b?,''cl_id'',?);',3:4)
c_eval('[Etot? eb_ang?]=irf_edb(e?(:,1:3),b?,90,''Epar'');',3:4)
c_eval('b_hat?=[diB?(:,1) diB?(:,2:4)./repmat(diB?(:,5),1,3)];',3:4)
c_eval('Epar?=irf_dot(Etot?,b_hat?);',3:4)

%c_eval('Epar?=[Epar?(:,1) (Epar?(:,2)-4)];',3:4)
c_eval('bxn?=irf_cross(b_hat?,[b_hat?(:,1) repmat([0 0 1],size(b_hat?,1),1)]);',3:4)
c_eval('Eper?=irf_dot(Etot?,bxn?);',3:4)

if 0
    c_eval('diPos?=c_coord_trans(''gse'',''dsi'',gsePos?,''cl_id'',?);',3:4)
    c_eval('diPospar?=irf_dot(diPos?,b_hat?);',3:4)
    dx_par=irf_add(1,diPospar3,-1,diPospar4);
    h_dx=sum(dx_par(2:(end-1),2))/(size(dx_par,1)-2);
end

[h_dt h_dt_error h_tplus h_tminus h_corr]=cn_mycorr(Epar3,Epar4,50,450,err);

vel=zdist/h_dt;

%cn_plot_corr(Epar3,Epar4,h_dt,'mycorr');
%%
vvec=[ones(size(Epar3,1),1) repmat(vel,size(Epar3,1),1)];
Tevec=[ones(size(Epar3,1),1) repmat(1/Teav,size(Epar3,1),1)];
c_eval('Phipar?t=irf_integrate(Epar?);',3:4)
c_eval('Phipar?=Phipar?t.*vvec;',3:4)

figure('name','Electron/ion hole');
np=4;
h=irf_plot(np);
isub=1;
if 1 % Electric field
    lb=[0.5 0.8 1]; % light blue
    lb=[0.6 0.7 1]; % vellight blue
    gr=[0 0.8 0]; % green
    
    c_eval('Epar?shift=Epar?; Epar?shift(:,1)=Epar?shift(:,1)+h_dt;',3:4);
    c_eval('Eper?shift=Eper?; Eper?shift(:,1)=Eper?shift(:,1)+h_dt;',3:4);
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,Epar3,'color',gr); hold(hca,'on');
    irf_plot(hca,Epar4,'b'); hold(hca,'on');
    irf_plot(hca,Epar4shift,'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,Epar3,'color',gr); hold(hca,'on');  
    irf_plot(hca,Epar4,'b'); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{||}  [mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,Eper3,'color',gr); hold(hca,'on');
    irf_plot(hca,Eper4,'b'); hold(hca,'on');
    irf_plot(hca,Eper4shift,'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,Eper3,'color',gr); hold(hca,'on');  
    irf_plot(hca,Eper4,'b'); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{\perp}  [mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]); 
end
if 1
    irf_filt()
end
if 1 % Potential and dB
    hca=h(isub); isub=isub+1;

    irf_plot(hca,Phipar3.*Tevec,'color',[0.8 0.5 0.2]); hold on;    
    irf_plot(hca,Phipar4.*Tevec,'color',[0.4 0.2 0.8]); hold on;
    irf_zoom('x',[cn_toepoch(t1) cn_toepoch(t2)]);
    set(hca,'ColorOrder',[[0.8 0.5 0.2];[0.4 0.2 0.8]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05])
    ylabel('phi_{||}/T_e')
end

%title(['Electrostatic potential and magnetic field fluctuations  \newline (n_e=',num2str(peaNeav,'%.3f'),' cm^{-3} B_0=',num2str(B0,'%.1f'),' nT)'])
irf_zoom(h,'x',[Epar3(1,1) Epar3(end,1)])

if 1 % Xtick and Xticklabel length
    xticks0=get(hca,'xtick');
    midindex=ceil(length(xticks0)/2);
    xticklabels=cell(1,length(xticks0));
    for l=1:length(xticks0); xticklabels{l}=' '; end 
    if mod(length(xticks0-1),4)==0; si=1; else si=2; end
    
    for l=si:2:length(xticks0)
        xticklabels{l}= num2str((xticks0(l)-xticks0(midindex))*vel,'%.0f');
    end
    set(h(1),'xaxislocation','top','color','white','xtick',xticks0,'xticklabel',xticklabels);    
    xlabel(h(1),'km');
end

if 1 
    Ld=irf_plasma_calc(B0,peaNeav,peaNeav,Teav,900,'Ld');
    delt=(Ld/1000)/abs(vel);
    dt1=Epar3(fix(end/2),1);
    dt2=dt1+delt*10;
    
    for k=1:np
        irf_pl_mark(h(k),[dt1 dt2],[1 1 0.0])
    end   
end


title(h(1),['\Delta t=',...
    num2str(h_dt,'%.4f'),' s   v=',...
    num2str(vel,'%.0f'),' km/s   T_e=',...
    num2str(Teav,'%.0f'),' eV   \lambda_{De}=',...
    num2str(Ld/1000,'%.1f'),' km  (yellow marking=10\lambda_{De})']);
%'Filtered and shited electric field, electrostatic potential and filtered magnetic field')
xy=get(gcf,'Position');
set(gcf,'Position',[xy(1) xy(2) xy(3)*1 xy(4)*1.2])

abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)'};
for k=1:np
    irf_legend(h(k),[abcde(k)],[0.02,0.92],'color',[0 0 0])
end



if 0 % Marking 10 lambda_De
    De=1;
    deltat=lambda/vchoice;
    tfirst=phi3(fix(end/2),1);

    tsecond=tfirst-deltat;
    for k=1:np
        irf_pl_mark(h(k),[tfirst tsecond],[1 1 0.0])
    end
    if 0
        tinterval=gsmB3(end,1)-gsmB3(1,1);
        ind=1;
        for k=1:size(gsmB3(end,1),1)
                vec(ind,1)=[gsmB3(k,1) 10]; 
                ind=ind+1;
          
        end
        irf_plot(h(2),vec)
    end
    %markvec(:,1)=[tfirst;tsecond];
    %markvec(:,2)=[60 60];
    %irf_plot(h(2),markvec,'+-')
end
set(gcf,'PaperPositionMode','auto');
%irf_zoom(h,'x',[Epar3(1,1) Epar3(end,1)])