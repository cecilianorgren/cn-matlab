

clear 
Tcomputsta=clock; 
%--------------------------------------------

Tsta='2003-10-02T00:46:54Z'; 
Tend='2003-10-02T00:47:07Z'; 
Tnull='2003-10-02T00:47:00.87Z'; 

%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
%--------------------------------------------

caa_load C


%background magnetic field
for ic=1:4
  c_eval(['B?=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');'],ic);
  c_eval(['R?=getmat(C?_CP_AUX_POSGSE_1M,''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'');'],ic);
  
  c_eval(['B?=irf_gse2gsm(B?);'],ic);
  c_eval(['R?=irf_gse2gsm(R?);'],ic);
end

B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);


for ic=1:4
  c_eval(['R?=irf_resamp(R?,B?);'],ic);
end
for ic=1:4
  c_eval(['B?=irf_tlim(B?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
end


%--------------------------------------------
%change the coordinates to eigen vector corrdinate
e1=[0.650, -0.727, 0.222]; e1=irf_norm(e1);
etmp=[0 1 0];
e3=cross(e1,etmp);          e3=irf_norm(e3);
e2=cross(e3,e1);

L=[0.79 0.58 0.19]; L=irf_norm(L);
etmp=[0 0 1];
M=cross(etmp,L);          M=irf_norm(M);
N=cross(L,M);

for ic=1:4,
    c_eval(['B?=irf_newxyz(B?, e1,e2,e3);'],ic);
    c_eval(['R?=irf_newxyz(R?, e1,e2,e3);'],ic);
end
%--------------------------------------------

gradB=c_4_grad('R?','B?','grad');
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');

%construct B around null
idxnull=find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
dB_null=reshape(gradB(idxnull,2:end),3,3);

     %[V,D] = eig(dB_null);
     %numda=[D(1,1) D(2,2) D(3,3)]
     %a1=['[' num2str(V(1,1),'%.3f') ', ' num2str(V(2,1),'%.3f') ', ' num2str(V(3,1),'%.3f') ']']
     %a2=['[' num2str(V(1,2),'%.3f') ', ' num2str(V(2,2),'%.3f') ', ' num2str(V(3,2),'%.3f') ']']
     %a3=['[' num2str(V(1,3),'%.3f') ', ' num2str(V(2,3),'%.3f') ', ' num2str(V(3,3),'%.3f') ']']

dR1=inv(dB_null) * (B1(idxnull,2:4))';
R_null=R1(idxnull,2:4)-dR1';

Rsc1=R1(idxnull,2:end);
Rsc2=R2(idxnull,2:end);
Rsc3=R3(idxnull,2:end);
Rsc4=R4(idxnull,2:end);

Rsc1=Rsc1-R_null;
Rsc2=Rsc2-R_null;
Rsc3=Rsc3-R_null;
Rsc4=Rsc4-R_null;

j_null=j(idxnull,2:end) .* 1e9;



%% construct B around null
set(0,'defaultLineLineWidth', 0.5);
set(0,'defaultAxesFontSize', 12);
set(0,'defaultTextFontSize', 12);
set(0,'defaultAxesFontUnits', 'pixels');
fig1=figure( ...
          'Name','Dataset coverage', ...
          'Tag','XYXgsm');clf;
set(fig1,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
xLeft = 0; yTop = -1;
set(fig1,'PaperPosition',[xLeft yTop xSize ySize])
set(fig1,'Position',[10 10 xSize*coef ySize*coef])


h=[];

h(1)=axes('position',[0.1 0.1 0.8 0.8]); % [x y dx dy]
aaa=1;

BoxWid=1500; 


for Xgrid=[-10 10]
    %for theta=[0:120:240]*pi/180
    for theta=[0:30:330]*pi/180
        Ygrid=10*cos(theta);
        Zgrid=10*sin(theta);
        
        %-----inverse trace-----
        ii=1;
        Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
        while abs(Xprev)<BoxWid & abs(Yprev)<BoxWid & abs(Zprev)<BoxWid
            Xcurt=Xprev;
            Ycurt=Yprev;
            Zcurt=Zprev;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
            step=Bmcurt;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xprev=Xcurt-stepvec(1);
            Yprev=Ycurt-stepvec(2);
            Zprev=Zcurt-stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
            plot3(gca, Xline, Yline, Zline); hold on;
            %Nlin=length(Xline);
            %cline(Xline, Yline, Zline, Bmline, 0, 40, jet); view(3); hold on;
            %arrP=fix(Nlin*0.9);
            %arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'b',2,5); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
        
        
        %-----positive trace-----
        ii=1;
        Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
        while abs(Xnext)<BoxWid & abs(Ynext)<BoxWid & abs(Znext)<BoxWid
            Xcurt=Xnext;
            Ycurt=Ynext;
            Zcurt=Znext;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
            step=Bmcurt;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xnext=Xcurt+stepvec(1);
            Ynext=Ycurt+stepvec(2);
            Znext=Zcurt+stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
            plot3(gca, Xline, Yline, Zline,'r'); hold on;
            %Nlin=length(Xline);
            %cline(Xline, Yline, Zline, Bmline, 0, 40, jet); view(3); hold on;
            %arrP=fix(Nlin*0.97);
            %arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'b',2,80); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
     
    end
end


% Jlinelngth=BoxWid;
% Jendp1=Jlinelngth*irf_norm(j_null);
% Jendp2=-Jendp1;
% arrow3(Jendp2, Jendp1,'s5',3,5); hold on;
% light('position',Jendp1); lighting gouraud;


plot3(gca, [Rsc1(1) Rsc2(1) Rsc3(1) Rsc4(1)], [Rsc1(2) Rsc2(2) Rsc3(2) Rsc4(2)], ...
           [Rsc1(3) Rsc2(3) Rsc3(3) Rsc4(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc2(1) Rsc4(1) Rsc1(1) Rsc3(1)], [Rsc2(2) Rsc4(2) Rsc1(2) Rsc3(2)], ...
           [Rsc2(3) Rsc4(3) Rsc1(3) Rsc3(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc1(1)], [Rsc1(2)],[Rsc1(3)], 'ks', 'Linewidth',1, ...
           'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12); hold on;
plot3(gca, [Rsc2(1)], [Rsc2(2)],[Rsc2(3)], 'rs', 'Linewidth',1, ...
           'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;
plot3(gca, [Rsc3(1)], [Rsc3(2)],[Rsc3(3)], 'gs', 'Linewidth',1, ...
           'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12); hold on;
plot3(gca, [Rsc4(1)], [Rsc4(2)],[Rsc4(3)], 'bs', 'Linewidth',1, ...
           'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12); hold off;   
 

% maxB=max(Bmax);
% minB=min(Bmin);
       
% hcb=colorbar('peer',h(1),'North', 'XDir','reverse', 'TickDir','out', 'XAxisLocation','top');
% posFig=get(h(1),'Position'); 
% left=posFig(1)+posFig(3)*0.7; low=posFig(2)+posFig(4)*4/5; width=posFig(3)/4; height=0.015;
% set(hcb,'Position',[left low width height]);
%ylabel(hcb,'|B|');


caxis([0, 40]);   %here is derived from minB & maxB. it should be consistent with cline range
set(gca,'Xlim',[-1500 1500], 'Ylim',[-1500 1500], 'Zlim',[-1500 1500]);
set(gca,'xtick',[-1500:500:1500], 'ytick',[-1500:500:1500], 'ztick',[-1500:500:1500]);
set(gca,'DataAspectRatio',[1 1.0 1]);
xlabel(gca,'e_{1} [km]');
ylabel(gca,'e_{2} [km]');
zlabel(gca,'e_{3} [km]');
grid off;



%angles=get(gca,'view');
%set(gca,'view',angles)



%% save figure
% set(fig1,'render','painters');
% figname=['B_around_null_cycle10000_spiral_cluster'];
% print(fig1, '-dpdf', [figname '.pdf']);


set(fig1,'renderer','opengl');
%axis off;
figname=['Btopology_3D'];
print(fig1, '-dpng','-r400',[figname '.png']);

set(gca,'view',[0 90])
figname=['Btopology_xy_plane'];
print(fig1, '-dpng','-r400',[figname '.png']);

set(gca,'view',[0 0])
figname=['Btopology_xz_plane'];
print(fig1, '-dpng','-r400',[figname '.png']);

set(gca,'view',[90 0])
figname=['Btopology_yz_plane'];
print(fig1, '-dpng','-r400',[figname '.png']);



Tcomputend=clock;
Telaps=etime(Tcomputend, Tcomputsta)/60  %minute