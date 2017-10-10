% Daniels reconnection events


events = '2';
switch events
    case '1' % 22 April 2008 - 17:00:00 UT - 18:30:00 UT
        tint = toepoch([2008 04 22 17 00 00;2008 04 22 18 30 00])';
        savePath = '/Users/Cecilia/Data/Cluster/20080422';
    case '2' % 17 April 2007 - 15:20:00 UT - 16:10:00 UT
        tint = toepoch([2007 04 17 15 20 00;2007 04 17 16 10 00])';
        savePath = '/Users/Cecilia/Data/Cluster/20070417';        
end
sclist = 1:4;
cd(savePath)
%%
for kk=7
    dataproducts = {'C?_CP_FGM_FULL_ISR2','C?_CP_EFW_L2_E3D_INERT','C?_CP_PEA_PITCH_SPIN_DEFLUX','C?_CP_STA_CWF_HBR_ISR2','C?_CP_STA_CWF_GSE','C?_CP_PEA_MOMENTS','C?_CP_CIS_HIA_ONBOARD_MOMENTS'};
    caa_download(tint,dataproducts{kk})
end