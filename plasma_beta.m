irf_units;
%% Ion particle pressure
[caaIP,~,codifIPH]=c_caa_var_get('pressure__C4_CP_CIS_CODIF_HS_H1_MOMENTS');
for k=1:size(codifIPH.data,1)
    IPH(k,1)=codifIPH.t(k);
    IPH(k,2)=sum(diag(squeeze(codifIPH.data(k,:,:,:))))/3;
end
%% Ion particle pressure
[caaIPO,~,codifIPO]=c_caa_var_get('pressure__C4_CP_CIS_CODIF_HS_O1_MOMENTS');
for k=1:size(codifIPO.data,1)
    IPO(k,1)=codifIPO.t(k);
    IPO(k,2)=sum(diag(squeeze(codifIPO.data(k,:,:,:))))/3;
end
%% Ion particle pressure
[caaIPHe,~,codifIPHe]=c_caa_var_get('pressure__C4_CP_CIS_CODIF_HS_He1_MOMENTS');
for k=1:size(codifIPHe.data,1)
    IPHe(k,1)=codifIPHe.t(k);
    IPHe(k,2)=sum(diag(squeeze(codifIPHe.data(k,:,:,:))))/3;
end
%% Ion particle pressure from N and T
[caaIN,~,codifINH]=c_caa_var_get('density__C4_CP_CIS_CODIF_HS_H1_MOMENTS');
[caaIT,~,codifITH]=c_caa_var_get('T__C4_CP_CIS_CODIF_HS_H1_MOMENTS');
IPNT=irf_multiply(Units.kB*1e6*1e6*1e9,codifINH,1,codifITH,1);
%% Ion particle pressure HIA
[caahiaIP,~,hiaIPH]=c_caa_var_get('pressure__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
%% Electron particle pressure
[caaEP,~,peaceEP]=c_caa_var_get('Data_Pressure_GSE__C4_CP_PEA_MOMENTS');
for k=1:size(peaceEP.data,1)
    EP(k,1)=peaceEP.t(k);
    EP(k,2)=sum(diag(squeeze(peaceEP.data(k,:,:,:))))/3;
end
%%
[caaIN,~,codifINO]=c_caa_var_get('density__C4_CP_CIS_CODIF_HS_O1_MOMENTS');
[caaEN,~,peaceEN]=c_caa_var_get('Data_Density__C4_CP_PEA_MOMENTS');
%% Magnetic field pressure
irf_units;
PB=[gsmB4fgm(:,1) gsmB4fgm(:,5).^2*1e-9/2/Units.mu0];

%%
Ptot=irf_add(1,irf_add(1,IPH,1,EP),1,PB);
%% Plasma beta
betaENT=irf_multiply(1,irf_add(1,IPNT,1,EP),1,PB,-1);
betaEH=irf_multiply(1,irf_add(1,IPH,1,EP),1,PB,-1);
betaH=irf_multiply(1,IPH,1,PB,-1);
%beta=irf_multiply(1,Ptot,1,PB,-1);
figure(29);irf_plot({betaENT,betaEH,betaH},'comp');set(gca,'yscale','log')
legend('e-, N and T (H+)','e- and H+','H+')
grid off;
%% Plot
figure(27);irf_plot({IPH,hiaIPH,IPHe,IPO,EP,PB,Ptot},'comp');
legend('H+ (CODIF)','H+ (HIA)',...
    'He+ (CODIF)','O+ (CODIF)',...
    'e- (PEACE)','Magnetic field (FGM)',...
    'H+,e-,BTotal pressure')
ylabel('[nPa]')
title('Pressure C4')
irf_zoom('x',[IPH(1,1) IPH(end,1)]);
set(gcf,'PaperPositionMode','auto');
%%
figure(71);irf_plot({[codifINH(:,1) codifINH(:,2)*1000],codifITH},'comp')