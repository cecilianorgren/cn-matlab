% Check for data and download
%caa_download(tailBMC3.TimeInterval(50,:),'list:C3*EFW*')
iEvent = 34; % Daniels event
tintDL = tailBMC4.TimeInterval(iEvent,:);
tintStr = irf_time(tintDL(1),'yyyymmdd');
savePath = ['/Users/Cecilia/Data/BM/' tintStr];
if ~isdir(savePath); mkdir(savePath); end
cd(savePath)
caa_download(tintDL,'C4_CP_EFW_L2_E3D_INERT')
caa_download(tintDL,'C4_CP_FGM_FULL_ISR2')
caa_download(tintDL,'C4_CP_PEA_PITCH_3DXH_PSD')
caa_download(tintDL,'C4_CP_PEA_3DXPH_PSD')
caa_download(tintDL,'C4_CP_PEA_3DXH_PSD')
caa_download(tintDL,'C4_CP_PEA_PITCH_3DXH_DEFlux')
