%% Plot fields
if 1 % Time for electron holes
    ta=[2007 08 31 10 10 1 1]; tb=[2007 08 31 10 20 1 1];
    
    
    ev='B2'; t1=[2007 08 31 10 17 46.45]; t2=[2007 08 31 10 17 46.70]; % the last from B
    
    ev='X'; t1=[2007 08 31 10 17 44.50]; t2=[2007 08 31 10 17 46.90];
    ev='X'; t1=[2007 08 31 10 17 51.60]; t2=[2007 08 31 10 17 51.90]; % one separate
    
    ev='X'; t1=[2007 08 31 10 18 51.10]; t2=[2007 08 31 10 18 51.70]; % dipuls
    ev='DL'; t1=[2007 08 31 10 18 42.50]; t2=[2007 08 31 10 18 43.00]; % dipuls
    ev='C1'; t1=[2007 08 31 10 17 54.10]; t2=[2007 08 31 10 17 54.60];
    ev='C2'; t1=[2007 08 31 10 17 52.60]; t2=[2007 08 31 10 17 52.95];
    ev='C4'; t1=[2007 08 31 10 17 51.85]; t2=[2007 08 31 10 17 52.05];
    ev='C3'; t1=[2007 08 31 10 17 51.65]; t2=[2007 08 31 10 17 51.85];
    ev='X'; t1=[2007 08 31 10 18 35.40]; t2=[2007 08 31 10 18 37.00];
    ev='D'; t1=[2007 08 31 10 17 36.00]; t2=[2007 08 31 10 17 36.60];
    ev='D'; t1=[2007 08 31 10 18 37.65]; t2=[2007 08 31 10 18 38.00];
    ev='A'; t1=[2007 08 31 10 18 35.80]; t2=[2007 08 31 10 18 36.50];
    ev='E1'; t1=[2007 09 02 15 21 20.50]; t2=[2007 09 02 15 21 21.50]; % many small pulses
    ev='E1'; t1=[2007 09 02 18 01 01.00]; t2=[2007 09 02 18 01 05.00]; % almost exact E3/E4 all the time
    ev='E1'; t1=[2007 09 02 18 01 11.00]; t2=[2007 09 02 18 01 13.00]; % almost exact E3/E4 all the time
    ev='F1'; t1=[2007  9 26 10 50 24.58];   t2=[2007  9 26 10 50.24 77];
    
    %ev='Y'; t1=[2007 08 31 10 16 03 00]; t2=[2007 08 31 10 16 05 00]; % seen only on one sc
    ev='B'; t1=[2007 08 31 10 17 45.50]; t2=[2007 08 31 10 17 46.70]; % many holes
    ev='B1'; t1=[2007 08 31 10 17 45.50]; t2=[2007 08 31 10 17 45.90]; % one hole
    ev='C'; t1=[2007 08 31 10 17 50.00]; t2=[2007 08 31 10 17 54.60]; % many holes
    ev='C2'; t1=[2007 08 31 10 17 50.00]; t2=[2007 08 31 10 17 53.50]; % many holes, half above
    ev='B2'; t1=[2007 08 31 10 17 46.45]; t2=[2007 08 31 10 17 46.65]; % one hole
    ev='B1'; t1=[2007 08 31 10 17 46.36]; t2=[2007 08 31 10 17 46.48]; % one hole
    ev='F1'; t1=[2007 08 31 10 18 07.0]; t2=[2007 08 31 10 18 07.40]; % one hole
end
if 1 % not eh's but waves
    ev='W1'; t1=[2007 08 31 10 17 20.50]; t2=[2007 08 31 10 17 21.50]; % 
    ev='W2'; t1=[2007 08 31 10 17 33.40]; t2=[2007 08 31 10 17 33.80]; % 
    ev='W3'; t1=[2007 08 31 10 17 36.00]; t2=[2007 08 31 10 17 36.60]; % 
    ev='W4'; t1=[2007 08 31 10 17 49.40]; t2=[2007 08 31 10 17 50.00]; % 
    ev='W5'; t1=[2007 08 31 10 18 07.80]; t2=[2007 08 31 10 18 08.60]; % 
    ev='W6'; t1=[2007 08 31 10 18 09.20]; t2=[2007 08 31 10 18 10.40]; % 
    ev='W7'; t1=[2007 08 31 10 18 16.20]; t2=[2007 08 31 10 18 17.20]; %
    ev='W8'; t1=[2007 08 31 10 18 18.80]; t2=[2007 08 31 10 18 19.30]; %
    ev='W9'; t1=[2007 08 31 10 18 19.80]; t2=[2007 08 31 10 18 20.60]; %
    ev='W10'; t1=[2007 08 31 10 18 35.80]; t2=[2007 08 31 10 18 36.50];
    ev='W11'; t1=[2007 08 31 10 18 37.30]; t2=[2007 08 31 10 18 37.70];
end
%%
tint= [cn_toepoch(t1) cn_toepoch(t2)];
if 1 % make date/time strings
    t1str_p=datestr(epoch2date(tint(1)),'yyyymmddTHHMMSSFFF');
    t2str_p=datestr(epoch2date(tint(2)),'MMSSFFF');
end
if 1 % Taking average density and temperature during time interval
    c_eval('peaNe?z=cn_toepoch(t1,t2,peaNe?);',3:4);
    c_eval('peaNe?av=sum(peaNe?z(:,2),1)/size(peaNe?z,1);',3:4);
    peaNeav=(peaNe3av+peaNe4av)/2;
    Ne=peaNeav*10^6; % m^-3
    
    c_eval('peaTepar?z=cn_toepoch(t1,t2,parTe?);',3:4);
    Teparav=mean([peaTepar3z(:,2);peaTepar4z(:,2)])*8.61734*10;
    c_eval('peaTeper?z=cn_toepoch(t1,t2,perTe?);',3:4);
    Teperav=mean([peaTeper3z(:,2);peaTeper4z(:,2)])*8.61734*10;
    Teav=(Teparav+Teperav)/2;
end

if 1 % Constants
    mu0=4*pi*1e-7;
    e=1.6e-19;
end
%
if 1 % Reconstructing Etot assuming it is parallel field
    c_eval('b?=cn_toepoch(t1,t2,gsmB?);',3:4)
    c_eval('diB?=c_coord_trans(''gsm'',''dsi'',b?,''cl_id'',?);',3:4)
    c_eval('[diEt? eb_ang?]=irf_edb(cn_toepoch(t1,t2,diE?),diB?,90,''Epar'');',3:4);
    %c_eval('diEt?=cn_toepoch(t1,t2,diE?);diEt?(:,4)=0;',3:4);
    c_eval('gsmEt?=c_coord_trans(''dsi'',''gsm'',diEt?,''cl_id'',?);',3:4)
    %EBang=mean([eb_ang3;eb_ang4]);
end
%
if 0 % Reconstructing Etot (longer interval) assuming it is parallel field
    c_eval('blong?=irf_tlim(gsmB?,cn_toepoch(ta),cn_toepoch(tb));',3:4)
    c_eval('diBlong?=c_coord_trans(''gsm'',''dsi'',blong?,''cl_id'',?);',3:4)
    c_eval('[diEtlong? eb_anglong?]=irf_edb(cn_toepoch(ta,tb,diE?),diBlong?,90,''Epar'');',3:4);
    c_eval('gsmEtlong?=c_coord_trans(''dsi'',''gsm'',diEtlong?,''cl_id'',?);',3:4)
    EBanglong=mean([eb_anglong3;eb_anglong4]);
end
%
if 1 % Transforming E and B to field aligned coordinates:
    % z along b, y 
    % perpendicular to b but in sc spin plane 
    % x third direction
    z=mean(irf_norm([b3(:,1:4); b4(:,1:4)]));
    %c_eval('zav?=[z?(:,1) repmat(sum(z?(:,2:4),1)/size(z?,1),size(z?,1),1)];',3:4);
    %c_eval('z?=zav?;',3:4);
    c_eval('sagsm?=c_coord_trans(''dsi'',''gsm'',[b?(:,1) repmat([0 0 1],size(b?,1),1)],''cl_id'',?);',3:4)
    sagsm=mean([sagsm3; sagsm4]);
    x=irf_norm(irf_cross(z,irf_cross(sagsm,z)));
    y=irf_norm(irf_cross(z,x));
    %c_eval('x?=irf_cross(z?,irf_cross(sagsm?,z?));',3:4)
    %c_eval('y?=irf_cross(z?,x?);',3:4)
    %xav=irf_norm(mean([x3(:,2:4); x4(:,2:4)],1));
    %yav=irf_norm(mean([y3(:,2:4); y4(:,2:4)],1));
    %zav=irf_norm(mean([z3(:,2:4); z4(:,2:4)],1));
    c_eval('facB?=[b?(:,1) b?(:,2:4)*x(2:4)'' b?(:,2:4)*y(2:4)'' b?(:,2:4)*z(2:4)''];',3:4)
    %c_eval('facE?=[gsmEt?(:,1) irf_dot(x,gsmEt?,1) irf_dot(y,gsmEt?,1) irf_dot(z,gsmEt?,1)];',3:4)
    c_eval('facE?=[gsmEt?(:,1) gsmEt?(:,2:4)*x(2:4)'' gsmEt?(:,2:4)*y(2:4)'' gsmEt?(:,2:4)*z(2:4)''];',3:4)
    %c_eval('facE?c=irf_lmn(gsmEt?,xav,yav,zav);',3:4)
end

if 0 % Transforming E and B to field aligned coordinates (longer interval):
    % z along b, y 
    % perpendicular to b but in sc spin plane 
    % x third direction
    c_eval('zlong?=irf_norm(blong?);',3:4);
    c_eval('zavlong?=[zlong?(fix(end/2),1) mean(zlong?(:,2:4),1)];',3:4);
    c_eval('zlong?=zavlong?;',3:4);
    c_eval('spgsmlong?=c_coord_trans(''dsi'',''gsm'',[blong?(:,1) repmat([1 1 0],size(blong?,1),1)/sqrt(2)],''cl_id'',?);',3:4)
    c_eval('ylong?=irf_cross(zlong?,irf_cross(spgsmlong?,zlong?));',3:4)
    c_eval('xlong?=irf_cross(ylong?,zlong?);',3:4)

    c_eval('facBlong?=[blong?(:,1) irf_dot(xlong?,blong?,1) irf_dot(ylong?,blong?,1) irf_dot(zlong?,blong?,1)];',3:4)
    c_eval('facElong?=[blong?(:,1) irf_dot(xlong?,gsmEtlong?,1) irf_dot(ylong?,gsmEtlong?,1) irf_dot(zlong?,gsmEtlong?,1)];',3:4)
end
if 1 % Take AC fields
    f_filt=1;
    c_eval('facBac?=irf_filt(facB?,f_filt,0,450,5);',3:4)
    c_eval('facEac?=irf_filt(facE?,f_filt,0,450,5);',3:4)
end

if 0 % Plot fields
    if 1 % Bac,Eac for C3/C4
        figure(6);h=irf_plot(4);
        irf_plot(h(1),facBac3);
        irf_legend(h(1),{'x','y: in spin plane','z: along B'},[0.02 0.1])
        ylabel(h(1),'B_{AC}   C3')
        irf_plot(h(2),facE3);
        irf_legend(h(2),{'x','y: in spin plane','z: along B'},[0.02 0.1])
        ylabel(h(2),'E   C3')
        irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)])
        irf_plot(h(3),facBac4);
        irf_legend(h(3),{'x','y: in spin plane','z: along B'},[0.02 0.1])
        ylabel(h(3),'B_{AC}   C4')
        irf_plot(h(4),facE4);
        irf_legend(h(4),{'x','y: in spin plane','z: along B'},[0.02 0.1])
        ylabel(h(4),'E   C4')
        irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)]);
        title(h(1),['B filtered at f_{filt}=',num2str(f_filt,'%.1f'),' Hz']);
        set(gcf,'PaperPositionMode','auto');
        eval(['print -dpng /Users/Cecilia/Dropbox/Cecilia/EH/20120820/',t1str_p,'-',t2str_p,'_Fields_',num2str(f_filt,'%.0f'),'Hz.png']);
    end
    if 1 % Like Tao C3
        figure(7);h=irf_plot(6);isub=1;
        hca=h(isub);irf_plot(hca,facE3(:,[1 4]));ylabel(hca,'E_{||}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac3(:,[1 4]));ylabel(hca,'\delta B_{||}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facE3(:,[1 2]));ylabel(hca,'E_{x}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac3(:,[1 3]));ylabel(hca,'\delta B_{y}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facE3(:,[1 3]));ylabel(hca,'E_{y}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac3(:,[1 2]));ylabel(hca,'h\delta B_{x}');isub=isub+1;
        irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)]);
        title(h(1),['All data C3, B filtered at f_{filt}=',num2str(f_filt,'%.1f'),' Hz']);
        set(gcf,'PaperPositionMode','auto');
        eval(['print -dpng /Users/Cecilia/Dropbox/Cecilia/EH/20120820/',t1str_p,'-',t2str_p,'_Fields3_',num2str(f_filt,'%.0f'),'Hz.png']);
    end
    if 1 % Like Tao C4
        figure(8);h=irf_plot(6);isub=1;
        hca=h(isub);irf_plot(hca,facE4(:,[1 4]));ylabel(hca,'E_{||}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac4(:,[1 4]));ylabel(hca,'\delta B_{||}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facE4(:,[1 2]));ylabel(hca,'E_{x}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac4(:,[1 3]));ylabel(hca,'\delta B_{y}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facE4(:,[1 3]));ylabel(hca,'E_{y}');isub=isub+1;
        hca=h(isub);irf_plot(hca,facBac4(:,[1 2]));ylabel(hca,'\delta B_{x}');isub=isub+1;
        irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)]);
        title(h(1),['All data C4, B filtered at f_{filt}=',num2str(f_filt,'%.1f'),' Hz']);
        set(gcf,'PaperPositionMode','auto');
        eval(['print -dpng /Users/Cecilia/Dropbox/Cecilia/EH/20120820/',t1str_p,'-',t2str_p,'_Fields4_',num2str(f_filt,'%.0f'),'Hz.png']); % 
    end
end

if 1 % Perform Time matching of parallel E, get dz
    % Get sc positions
    c_eval('gsmPos?=c_coord_trans(''gse'',''gsm'',cn_toepoch(t1,t2,gsePos?),''cl_id'',?);',3:4)
    dz=irf_dot(irf_add(1,gsmPos4,-1,gsmPos3),z);
    dzav=sum(dz(:,2)/size(dz,1));
    dy=irf_dot(irf_add(1,gsmPos4,-1,gsmPos3),y);
    dx=irf_dot(irf_add(1,gsmPos4,-1,gsmPos3),x);
    dyav=mean(dy(:,2));
    dxav=mean(dx(:,2));
    facdx=[dxav dyav dzav]*x(2:4)';
    facdy=[dxav dyav dzav]*y(2:4)';
    facdz=[dxav dyav dzav]*z(2:4)';
end
%
if 1 % Do time matching, get v
    facEac3=irf_resamp(facEac3,facEac4);
    window=50;err=0.05;
    %[dtc dt_error tplus tminus corr]=cn_mycorr(facEac3(:,[1 4]),facEac4(:,[1 4]),window,450,err);
    [dtx corrx]=cn_xcorr(facEac3,facEac4,window,'z');
    dt=-dtx;
    v=dz(:,2)/dt;
    vav=sum(v)/size(v,1);
end
if 1 % Integrate potential and do shift
    c_eval('Phi?=irf_integrate(facE?(:,[1 4]));',3:4)
    c_eval('Phi?(:,2)=Phi?(:,2)*vav;',3:4)
    Phishift4=Phi4; Phishift4(:,1)=Phishift4(:,1)-dt;
    Eshift4=facE4; Eshift4(:,1)=Eshift4(:,1)-dt;
end
if 1 % Calculate parallel current 
    %Niav=cn_toepoch(t1,Ni3);
    %Niav=Niav*1e6;
    N=peaNe4av*1e6;
    deltaB=irf_add(1,facB4,-1,facB3);
    j=cn_j(deltaB,[facdx facdy facdx]); % in A
    j(:,2:4)=j(:,2:4)*1e9;  % in nA
    %jz=irf_add(1,irf_multiply(10^(-3)/mu0,deltaB(:,[1 3]),1,facdx,-1),...
    %    -1,irf_multiply(10^(-3)/mu0,deltaB(:,[1 2]),1,[facdy],-1)); % nA
    evel=[j(:,1) j(:,4)*10^(-9)/N/e]; % km/s
    figure(10);h=irf_plot(3);isub=1;
    hca=h(isub);isub=isub+1;irf_plot(hca,deltaB);ylabel(hca,'\Delta B  [nT]'):legend(hca,'x','y','z')
    hca=h(isub);isub=isub+1;irf_plot(hca,j);ylabel(hca,'j  [nA/m^2]'):legend(hca,'x','y','z')
    %hca=h(isub);isub=isub+1;irf_plot(hca,irf_filt(jz,f_filt,0,450,5));ylabel(hca,['j_{||,AC} ',num2str(f_filt,'%.0f'),'  [nA]'])
    hca=h(isub);isub=isub+1;irf_plot(hca,evel);ylabel(hca,'(v_{i}-v_e)  [m/s]')
    title(h(1),['dx = ',num2str(facdx,'%.0f'),'  dy = ',num2str(facdy,'%.0f'),'  dz = ',num2str(facdz,'%.0f'),'         T_e=___eV, v_{te}=_____30*10^6 m/s, n=0.04 cc'])
    eval(['print -dpng 20120828/',t1str_p,'-',t2str_p,'_j.png']); %/Users/Cecilia/Dropbox/Cecilia/EH/20120820
end
if 1 % Plot shift and potential
    
    np=3;
    figure(9); set(gcf,'Position',[10 10 500 320]);h=irf_plot(np);isub=1;
    height=0.25; width=0.70;
    set(h(3),'Position',[0.17 0.1 width height]);
    set(h(2),'Position',[0.17 0.1+height width height])
    set(h(1),'Position',[0.17 0.1+2*height width height])
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facE3(:,[1 4]),facE4(:,[1 4]),facE4(:,[1 4])},...
        'comp','dt',[0 0 dt],'linestyle',{'-g','-b','--b'})   
    ylabel(hca,'E_{||} [mV/m]'); 
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 1]]);
    irf_legend(hca,{'C3','C4','C4_{shift}--'},[0.02 0.06])

    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facE3(:,[1 3]),facE4(:,[1 3]),facE4(:,[1 3])},...
        'comp','dt',[0 0 dt],'linestyle',{'-g','-b','--b'}) 
    ylabel(hca,'E_{y} [mV/m]'); 
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 1]]);
    irf_legend(hca,{'C3','C4','C4_{shift}--'},[0.02 0.06])
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{Phi3(:,[1 2]),Phi4(:,[1 2]),Phi4(:,[1 2])},...
        'comp','dt',[0 0 dt],'linestyle',{'-g','-b','--b'})  
    ylabel(hca,'\phi_{||} [V]'); 
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 1]]);
    irf_legend(hca,{'C3','C4','C4_{shift}--'},[0.02 0.06])
    %set(h(3),'ylim',[-200,700])
    
    if 1 % add second yaxis with phi/Te
        ax1 = gca;
        yticks_phi=get(ax1,'ytick');
        yticks_phi_new=yticks_phi/Teparav;
        ax2 = axes('Position',get(ax1,'Position'));
        set(ax2,'XAxisLocation','top',...
            'YAxisLocation','right','box','off',...
            'ytick', get(ax1,'ytick'),...
            'yticklabel',num2str(get(ax1,'ytick')'/Teparav,'%.2f'),...
            'ylimmode','manual','ylim',get(ax1,'ylim'),...
            'xtick',[],'Color','none');
        ylabel(ax2,'\phi/T_e');
    end
    
    if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,Phi3(:,[1 2]),'g'); hold(hca,'on');
    irf_plot(hca,Phi4(:,[1 2]),'b'); hold(hca,'on');
    irf_plot(hca,Phishift4(:,[1 2]),'b--'); hold(hca,'on');
    ylabel(hca,'\phi_{||} [V]'); 
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 1]]);
    irf_legend(hca,{'C3','C4','C4_{shift}--'},[0.02 0.04])
    end
    
    %hca=h(isub);isub=isub+1;
    %irf_plot(hca,j(:,[1 4]));ylabel(hca,'j_{||}  [nA/m^2]')
    
    %hca=h(isub);isub=isub+1;irf_plot(hca,evel);ylabel(hca,'(v_{i}-v_e)_{||}  [m/s]')
    
    
    irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)])

    if 1 % Mark 10 Debye lengths
        B0=sum(b3(:,5)/size(b3,1));
        Ld=irf_plasma_calc(B0,peaNeav,peaNeav,Teav,900,'Ld');
        delt=(Ld/1000)/abs(sum(v,1)/size(v,1));
        dt1=facEac3(find(facEac3(:,4)==max(facEac3(:,4))),1)-0.014;
        
        dt2=dt1+delt*10;   
        for k=1:np
            irf_pl_mark(h,[dt1 dt2],'yellow')
        end   
    end
    for k=1:3; grid(h(k),'off'); end

    title(h(1),['\Delta t=',...
        num2str(dt*1000,'%.1f'),' ms   v=',...
        num2str(vav,'%.0f'),' km/s   T_{e,||}=',...
        num2str(Teparav,'%.0f'),' eV   \lambda_{De}=',...
        num2str(Ld/1000,'%.1f'),' km  (yel mark=10\lambda_{De})']);

    if 1 % Add parellel length scale on top
        xticks0=get(hca,'xtick');
        midindex=ceil(length(xticks0)/2);
        xticklabels=cell(1,length(xticks0));
        for l=1:length(xticks0); xticklabels{l}=' '; end 
        if mod(length(xticks0-1),4)==0; si=1; else si=2; end

        for l=si:2:length(xticks0)
            xticklabels{l}= num2str((xticks0(l)-xticks0(midindex))*vav,'%.0f');
        end
        set(h(1),'xaxislocation','top','color','white','xtick',xticks0,'xticklabel',xticklabels);    
        xlabel(h(1),'km');
    end   
    %irf_pl_mark(h,[dt1 dt2],'yellow')
    xy=get(gcf,'Position');
    set(gcf,'Position',[xy(1) xy(2) xy(3)*1 xy(4)*1.5])    
    set(gcf,'PaperPositionMode','auto');
    %irf_pl_mark(h,[dt1 dt2],'yellow')
    eval(['print -dpng /Users/Cecilia/Konferenser&Skolor/BoulderPostCluster12/',t1str_p,'-',t2str_p,'_Phi_E.png']); %/Users/Cecilia/Dropbox/Cecilia/EH/20120820
end

%% Some additional data that might be interesting
 Ti=mean(eVTi3(:,2));
 Bav=mean(diB3(:,5));
 Vith=sqrt(192*Ti);
 Veth=600*sqrt(Teav);
 CS=Vith*sqrt(Teav/Ti);
 
%% Create structure to save variables in
variables.x_hat=xav;
variables.y_hat=yav;
variables.z_hat=zav;
variables.dx=dxav;
variables.dy=dyav;
variables.dz=dzav;
variables.dt=dt;
variables.v=vav;
variables.T_i=Ti;
variables.T_e=Teav;
variables.B=Bav;
variables.V_ith=Vith;
variables.V_eth=Veth;
variables.C_s=CS;
variables.L_Debye=Ld/1000;
variables

%% Plot particles
if 0 % Par and perp electric field spectrograms summed over about a second..
    fs=450; overlap=20; nfft=512;
    c_eval('facEz?fft=irf_powerfft(cn_toepoch(t1,t2,facE?(:,[1 4])),nfft,fs,overlap);',3:4);
    c_eval('facEy?fft=irf_powerfft(cn_toepoch(t1,t2,facE?(:,[1 3])),nfft,fs,overlap);',3:4);
    
    facEz3fft.f_label=['f [Hz] z'];
    facEy3fft.f_label=['f [Hz] y'];
    figure(12);h=irf_plot(2);isub=1;
    hca=h(isub);isub=isub+1;irf_spectrogram(hca,facEz3fft);
    hca=h(isub);isub=isub+1;irf_spectrogram(hca,facEy3fft);
end

if 1 % Look at particle data, pitch angle distribution etc
    isub=1; ic=3; 
    if ic==3; figure(12); elseif ic==4; figure(15); end 
    h=irf_plot(7); 
    if 1 % Plot electric field
        hca=h(isub);isub=isub+1;
        c_eval('irf_plot(hca,facE?);',ic);
        irf_legend(hca,{'x','y: in spin plane','||'},[0.02 0.04])
    end
    if 0 % PEACE energy spectrogram (with irf_spectrogram)
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
        [var,dobj,varmat,varunits]=c_caa_var_get(varname);
        energies=varmat.dep_x{2}.data; % 2166x44 energy
        phi=varmat.dep_x{1}.data; % 2166x12 pitch angles
        flux=varmat.data; % particle flux
        time=varmat.t;
        flux(isnan(flux))=0.000;
        fluxsum0=sum(flux,2);
        fluxsum=zeros(size(fluxsum0,1),size(fluxsum0,3));

        for k=1:44;
            fluxsum(:,k)=fluxsum0(:,1,k);
        end
        fluxsum(:,26:end)=[];
        energies(:,26:end)=[];
        energies(isnan(energies))=0;
        eng=log10(fluxsum)+log10(energies);
        eng2=10.^eng;
        eng2(isinf(eng2))=9;

        %h=irf_plot(2);
        pc=pcolor(time(:,1),energies(1,:),eng'); shading flat; 
        %hcb = colorbar('peer',pc);
        %ylabel(hcb,{'log_{10} dEF';'keV/cm^2s sr keV'},'fontsize',12)

        irf_colormap(gca,'space');
        irf_timeaxis
        set(gca,'yscale','log','TickDir','Out')

        %%%%%%%%%%%%%%%%%

        if 0
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,{['C' num2str(3)]},[0.98 0.05],'color','k');
        ylabel('E [eV]');
        end
    end
    if 0 % PEACE energy spectrogram (with pcolor)
        hca=h(isub); isub=isub+1;
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
        [var,dobj,varmat,varunits]=c_caa_var_get(varname);
        energies=varmat.dep_x{2}.data; % 2166x44 energy
        phi=varmat.dep_x{1}.data; % 2166x12 pitch angles
        flux=varmat.data; % particle flux
        time=varmat.t;
        flux(isnan(flux))=0.000;
        fluxsum0=sum(flux,2);
        fluxsum=zeros(size(fluxsum0,1),size(fluxsum0,3));

        for k=1:44;
            fluxsum(:,k)=fluxsum0(:,1,k);
        end
        fluxsum(:,26:end)=[];
        energies(:,26:end)=[];
        energies(isnan(energies))=0;
        eng=log10(fluxsum)+log10(energies);
        eng2=10.^eng;
        eng2(isinf(eng2))=9;

        %h=irf_plot(2);
        ud = get(gcf,'userdata');
        t_st_e = double(ud.t_start_epoch);
        pc=pcolor(hca,double(time(:,1)-t_st_e),double(energies(1,:)),double(eng')); shading flat; 
        shading(hca,'flat')
        %hcb = colorbar('peer',pc);
        %ylabel(hcb,{'log_{10} dEF';'keV/cm^2s sr keV'},'fontsize',12)
        irf_legend(hca,{[varname,'   C' num2str(ic)]},[0.98 0.85],'color','k')
        irf_colormap(gca,'space');
        irf_timeaxis
        set(gca,'yscale','log','TickDir','Out')

        %%%%%%%%%%%%%%%%%

        if 0
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,{['C' num2str(3)]},[0.98 0.05],'color','k');
        ylabel('E [eV]');
        end
    end
    if 1 % RAPID electron spectrogram
        hca=h(isub); isub=isub+1;
        c_eval('irf_plot(hca,''Electron_Dif_flux__C?_CP_RAP_ESPCT6'',''colorbarlabel'',''log10 dF\newline 1/cm^2 s sr keV'',''fitcolorbarlabel'');',ic);
        caxis([0.51 4.49]);
        irf_legend(hca,{['C' num2str(ic)]},[0.98 0.85],'color','k')
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
    end
    if 1 % PEACE 3DXPH_DEFlux high res energy spectrogram
        hca=h(isub);isub=isub+1;
        res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_PSD',ic));
        [~,ind]=irf_tlim(res.tt,tint);
        specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
        if 1, % energy spectorgram (integrated over pitch angles)
            specrec.f=log10(res.en);
            specrec.p=res.omni(ind,:);
            specrec.f_label=['Log10 ' res.enlabel];
            irf_spectrogram(hca,specrec);
        elseif 0, % pitch angle spectrogram for given energy
            specrec.f=res.theta;specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            enindex=13;
            specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
            specrec.p=log10(res.data(ind,:,enindex));
            irf_spectrogram(hca,specrec);
            set(hca,'ytick',[30 60 90 120 150]);
        end
        caxis(hca,[-1.99 0.49]);
        irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
    end
    if 1 % Plot pitch angle data
        hca=h(isub); isub=isub+1;
        
        if ic==3; cor=4; elseif ic==4; cor=0; end
        
        res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',ic));
        [~,ind]=irf_tlim(res.tt,tint);
        specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        specrec.en=res.en(:);
        flabel=specrec.f_label;
        
        enindex=11;
        specrec.f_label=[flabel '   C' num2str(ic)   '\newline[E=' num2str(res.en(enindex-cor),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        irf_spectrogram(hca,specrec);
        irf_legend(hca,['C' num2str(ic)],[0.98 0.98])
        
        enindex=13;
        specrec.f_label=[flabel '   C' num2str(ic)   '\newline[E=' num2str(res.en(enindex-cor),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        hca=h(isub);isub=isub+1;
        irf_spectrogram(hca,specrec);
        irf_legend(hca,['C' num2str(ic)],[0.98 0.98])
        
        enindex=15;
        specrec.f_label=[flabel    '\newline[E=' num2str(res.en(enindex-cor),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        hca=h(isub);isub=isub+1;
        irf_spectrogram(hca,specrec);
        irf_legend(hca,['C' num2str(ic)],[0.98 0.98])
    end
    if 1 % Plot velocity
        c_eval('[caaparVe?,~,gseparVe?]=c_caa_var_get(''Data_Velocity_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'');',3:4);
        c_eval('[caaperVe?,~,gseperVe?]=c_caa_var_get(''Data_Velocity_ComponentPerpendicularToMagField__C?_CP_PEA_MOM'');',3:4);
        %c_eval('gseparVe?=irf_abs(gseparVe?);',3:4)
        
        c_eval('gsmparVe?=c_coord_trans(''gse'',''gsm'',gseparVe?,''cl_id'',?);',3:4)
        c_eval('parVe?=irf_dot(z?,gsmparVe?);',3:4)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,{parVe3,parVe4},'comp')
        irf_legend(hca,{'C3','C4'},[0.98 0.98])
    end
    irf_plot_axis_align
    irf_zoom(h,'x',[cn_toepoch(ta) cn_toepoch(tb)])
end

%% Plot pitch angle dta for many different energies

for ic=3:4
    %ic=4;
    if ic==3; cor=4; cor=0; cor2=0; fign=16; elseif ic==4; cor=0; cor2=4; fign=17; end
    
    % Construct data
    res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',ic));
    [~,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    specrec.en=res.en(:);
    flabel=specrec.f_label;
    
    
    n_pl=12;
    figure(fign);h=irf_plot(n_pl);isub=1;
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,facE?);',ic)
    step=2;
    start_index=5;
    for k=(start_index:step:(start_index+step*(n_pl-2)))+cor2
        hca=h(isub);isub=isub+1;
        specrec.f_label=[flabel '   C' num2str(ic)   '\newline[E=' num2str(res.en(k-cor),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,k));
        irf_spectrogram(hca,specrec);
        irf_legend(hca,['C' num2str(ic)],[0.98 0.98])
        caxis(hca,[4.5 7.5])
    end
    irf_plot_axis_align
    irf_zoom(h,'x',[facE3(1,1) facE3(end,1)])
    set(gcf,'PaperPositionMode','auto');
end

