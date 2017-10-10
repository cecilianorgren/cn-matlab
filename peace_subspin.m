




%Notice: load the CAA data first (CAA_LOAD) 
%--------------------------------------------

ic=2;
Tsta='2003-09-19T23:31:10.000Z';   
Tend='2003-09-19T23:31:50.000Z';

%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
yr=substr(Tsta,1,4); mo=substr(Tsta,6,2); dy=substr(Tsta,9,2);
hr_sta=substr(Tsta,12,2); mi_sta=substr(Tsta,15,2); sc_sta=substr(Tsta,18,2); 
hr_end=substr(Tend,12,2); mi_end=substr(Tend,15,2); sc_end=substr(Tend,18,2); 
date=[yr mo dy]; 
%--------------------------------------------


res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DRH_PSD',ic));
[delmett,ind]=irf_tlim(res.tt,tint);
specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  



%% Init figure
n_subplots=5;

i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(61);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = 1; yTop = -1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


lnwid=1; %0.75;


%% Bz plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
specrec.f=log10(res.en);
specrec.p=res.omni(ind,:);
specrec.f_label=['Log10 ' res.enlabel];
    
irf_spectrogram(gca,specrec);


% Bx plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
 specrec.f=res.theta;specrec.f_label='Pitch angle';
 specrec.p=res.pitch_angle(ind,:);
 enindex=15;
 res.en(enindex)
 specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
 specrec.p=log10(res.data(ind,:,enindex));

 irf_spectrogram(gca,specrec);

 % Bx plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
 enindex=14;
 res.en(enindex)
 specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
 specrec.p=log10(res.data(ind,:,enindex));

 irf_spectrogram(gca,specrec);
 

 
% Annotation
irf_zoom(tint,'x',h);
set(h(1:i_subplot-2),'XTickLabe','');

leg = 'abcdefghijklmn';
for ii=1:i_subplot-1
    set(h(ii),'Ylimmode','auto')
    set(h(ii),'ColorOrder',[0 0 0]);
    irf_legend(h(ii),{[leg(ii)]}, [0.02, 0.05], 'FontSize',20);    % [1.02, 0.7]
end

    

%% save
set(gcf,'render','painters');
figname=['PEACE_subspin'];
print(gcf, '-dpdf', [figname '.pdf']);









