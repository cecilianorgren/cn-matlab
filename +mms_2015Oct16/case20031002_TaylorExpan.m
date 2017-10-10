
thresold=0.25;    %for indentify X type or O type null
%the interval should not too long, otherwise the figure is not readable
%--------------------------------------------

Tsta='2003-10-02T00:46:00Z'; 
Tend='2003-10-02T00:48:00Z'; 

Tsta='2015-10-16T10:33:20.00Z'; 
Tend='2015-10-16T10:34:00.00Z';

Tsta='2015-10-16T10:33:45.30Z'; 
Tend='2015-10-16T10:33:48.00Z';

%tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); 

%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
%--------------------------------------------

% caa_download(tint,'C?_CP_FGM_FULL');
% caa_download(tint,'C?_CP_AUX_POSGSE_1M');

% caa_load C


%background magnetic field
sclist = 1:4;
c_eval('B?=dmpaB?brst;',sclist);
c_eval('R? = gseR?.resample(dmpaB?brst);',sclist)
%c_eval(['R?=getmat(C?_CP_AUX_POSGSE_1M,''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'');'],sclist);
%c_eval(['B?=irf_gse2gsm(B?);'],sclist);
%c_eval(['R?=irf_gse2gsm(R?);'],sclist);

c_eval('B? = [double(B?.time.epochUnix) double(B?.data)];',sclist)
c_eval('R? = [double(R?.time.epochUnix) double(R?.data)];',sclist)


indices=c_fgm_poincare_index(B1,B2,B3,B4);


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

for ic=1:4
  c_eval(['BmagC?=irf_abs(B?);'],ic);
end



gradB=c_4_grad('R?','B?','grad');

d_dot_B=c_4_grad('R?','B?','div');
d_cros_B=c_4_grad('R?','B?','curl');


%error of curolmeter
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

temp=irf_abs(d_cros_B);
d_cros_B_mag=temp(:,[1 5]);
err_curlmeter=irf_multiply(1,d_dot_B,1,d_cros_B_mag,-1);
err_curlmeter(:,2)=abs(err_curlmeter(:,2))*100;


%Null type identification
for ii=1:length(d_dot_B(:,1))
    mksizSim(ii)=4;
    
    deltB_null=reshape(gradB(ii,2:end),3,3);
    [V,D] = eig(deltB_null);
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='w';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            type(ii)='X'; clr(ii)='k'; faceclr(ii)='w';
        end
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='r';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            type(ii)='o'; clr(ii)='k'; faceclr(ii)='w';
        end
    end
    %=========================================================
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='w';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if min(abs([D(1,1) D(2,2) D(3,3)]))/max(abs([D(1,1) D(2,2) D(3,3)]))<thresold
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
                typeSim(ii)='X'; clrSim(ii)='b'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='X'; clrSim(ii)='r'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            typeSim(ii)='X'; clrSim(ii)='k'; faceclrSim(ii)='w'; mksizSim(ii)=7;
        end
        %------------------------------------------------------------------
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='r';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))/max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) < (thresold*2)
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
                typeSim(ii)='o'; clrSim(ii)='b'; faceclrSim(ii)='w';
            end
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='o'; clrSim(ii)='r'; faceclrSim(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            typeSim(ii)='o'; clrSim(ii)='k'; faceclrSim(ii)='w';
        end
        %------------------------------------------------------------------
    end
    %=========================================================
    
    eigVal_err(ii,2)=abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
    sumeigVal(ii,2)=D(1,1)+D(2,2)+D(3,3);
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;
end
eigVal_err(:,1)=d_dot_B(:,1);
sumeigVal(:,1)=d_dot_B(:,1);
eigVal_err_v2(:,1)=d_dot_B(:,1);


%find null position
for ii=1:length(B1(:,1))
    dBeach=reshape(gradB(ii,2:end),3,3);
    dR1(ii,2:4)=B1(ii,2:4)*inv(dBeach');
    dR2(ii,2:4)=B2(ii,2:4)*inv(dBeach');
    dR3(ii,2:4)=B3(ii,2:4)*inv(dBeach');
    dR4(ii,2:4)=B4(ii,2:4)*inv(dBeach');
end
dR1(:,1)=B1(:,1);
dR2(:,1)=B1(:,1);
dR3(:,1)=B1(:,1);
dR4(:,1)=B1(:,1);

dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1)=dRmag1(:,1);


% Tnull='2003-08-17T16:41:55.65Z'; 
% idxnull=find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
% deltB_null=reshape(gradB(idxnull,2:end),3,3)
% d_dot_B_null=d_dot_B(idxnull,2)
% d_cros_B_null=d_cros_B(idxnull,2:end);
% uncertanty=err_4C(idxnull,2)
% dist=dRmin(idxnull,2)
% eigVal_uncertain=eigVal_err(idxnull,2)
% [V,D] = eig(deltB_null);
% numda=[D(1,1) D(2,2) D(3,3)]
% a1=['[' num2str(V(1,1),'%.3f') ', ' num2str(V(2,1),'%.3f') ', ' num2str(V(3,1),'%.3f') ']']
% a2=['[' num2str(V(1,2),'%.3f') ', ' num2str(V(2,2),'%.3f') ', ' num2str(V(3,2),'%.3f') ']']
% a3=['[' num2str(V(1,3),'%.3f') ', ' num2str(V(2,3),'%.3f') ', ' num2str(V(3,3),'%.3f') ']']



pasue=1;


%% Init figure
n_subplots=8;

i_subplot=1;
set(0,'DefaultAxesFontSize',9);
set(0,'DefaultLineLineWidth', 1);
fn=figure(61);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = 5; yTop = -1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

Re=6372;
lnwid=0.5;


%% Bz plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(BmagC1(:,[1 5]), 'color','k', 'Linewidth',lnwid); hold on;
irf_plot(BmagC2(:,[1 5]), 'color','r', 'Linewidth',lnwid); hold on;
irf_plot(BmagC3(:,[1 5]), 'color','g', 'Linewidth',lnwid); hold on;
irf_plot(BmagC4(:,[1 5]), 'color','b', 'Linewidth',lnwid); hold off;

ylabel('Bz [nT]');
set(gca,'Ylim',[-2 20]) 
%grid off;


%% distance plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(dRmag1(:,[1 5]), 'color','k', 'Linewidth',lnwid); hold on;
irf_plot(dRmag2(:,[1 5]), 'color','r', 'Linewidth',lnwid); hold on;
irf_plot(dRmag3(:,[1 5]), 'color','g', 'Linewidth',lnwid); hold on;
irf_plot(dRmag4(:,[1 5]), 'color','b', 'Linewidth',lnwid); hold off;

set(gca, 'Ylim',[0 200]); % ,'Ytick',[0:500:2000]
ylabel('|r| [km]');
grid off;


%% minmium distance plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(gca, dRmin, 'color','k', 'Linewidth',lnwid); hold on;
for ii=1:length(dRmin(:,1))
    irf_plot(gca, dRmin(ii,:), [typeSim(ii) clrSim(ii)], 'MarkerSize',mksizSim(ii),'MarkerFaceColor',faceclrSim(ii), 'Linewidth',lnwid); hold on;
end
hold off

set(gca, 'Ylim',[0 200]); % ,'Ytick',[0:500:2000]
ylabel('|r|_{min} [km]');
grid off;


%% error plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([err_4C(:,1) err_4C(:,2)],'k', 'Linewidth',lnwid); hold on;
irf_plot([err_4C(:,1) err_4C(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;

grid off;
set(gca,'Ylim',[0 100]);
ylabel('|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}| [%]');


%% error eigen value
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
%irf_plot([eigVal_err(:,1) eigVal_err(:,2)], 'color','b', 'Linewidth',lnwid); hold on;
irf_plot(eigVal_err_v2, 'b', 'Linewidth',0.5); hold on;
irf_plot([eigVal_err(:,1) eigVal_err(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;

grid off;
set(gca,'Ylim',[0 100]);
ylabel('|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]');



%% Annotation
irf_adjust_panel_position
irf_zoom(tint,'x',h)

pause=1;

set(gcf,'render','painters');
date=[Tsta(1:4) Tsta(6:7) Tsta(9:10) '-' Tsta(12:13) Tsta(15:16) Tsta(18:19) '_' Tend(12:13) Tend(15:16) Tend(18:19)];
figname=['case20031002_TaylorExpan'];
print(gcf, '-dpdf', [figname '.pdf']);



