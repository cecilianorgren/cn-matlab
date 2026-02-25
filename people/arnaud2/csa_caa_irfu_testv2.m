function csa_caa_irfu_test(test_location)
if nargin==0, test_location='IRF-Uppsala-Sweden';end
top_message='Now performing IRF-U test from CAA & CSA  1.1.1. main server (release May 6 2014)';
disp(top_message);
log_filename=['IRF-U_test_from_',num2str(test_location),'_',num2str(date),'.txt'];
fid=fopen(num2str(log_filename),'w');
line1=['Dataset\t\t\tStart time\t\tEnd time\t\tCAA Bytes \tCSA Bytes\tCAA time\tCSA time\tt_CSA/t_CAA\tFile size difference (percentage)\n'];
line12save='Dataset\t\t\t\tStart time\t\tEnd time\t\tCAA Bytes \tCSA Bytes\tCAA time\tCSA time\tt_CSA/t_CAA\tFile size difference (percentage)';
line2='***********			***********		***********		***********	***********	***********	***********	***********	***********\n';
fprintf(fid,'%s \r\n',sprintf(line1(1:end-2)));
fprintf(fid,'%s \r\n',sprintf(line2(1:end-2)));

cprintf(line1);
cprintf(line2);
%C3_CP_EFW_L3_E3D_INERT                            630369          630288          7.699373        5.914979        0.768242        
%csa_caa_irfu_test_core('EFW',test_location,3,1,'2005-01-29T17:23:00Z','2005-01-30T08:23:00Z',fid);
 [temp1,temp2]=csa_caa_irfu_test_core('EFW',test_location,3,1,'2005-01-29T10:32:00Z','2005-01-30T08:32:00Z',fid);
file_size(1)=temp1;time_ratio(1)=temp2;

%C3_CP_PEA_3DXPH_cnts                    2005-01-29T04:42:00Z     2005-01-29T04:42:00Z     194082          194932          7.441671        3.265106        0.438760        
%csa_caa_irfu_test_core('PEACE',test_location,3,1,'2005-01-29T04:42:00Z','2005-01-29T14:42:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('PEACE',test_location,3,1,'2005-01-29T10:58:00Z','2005-01-29T17:58:00Z',fid);
file_size(2)=temp1;time_ratio(2)=temp2;

%C3_CP_RAP_PAD_L3DD                      2005-01-29T06:02:00Z     2005-01-29T08:02:00Z     353166          353120          7.333393        4.820794        0.657376        
%csa_caa_irfu_test_core('RAPID',test_location,3,1,'2005-01-29T06:02:00Z','2005-01-29T08:02:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('RAPID',test_location,3,1,'2005-01-29T16:39:00Z','2005-01-30T04:39:00Z',fid);
file_size(3)=temp1;time_ratio(3)=temp2;

%C3_CP_WHI_ELECTRON_DENSITY              2005-01-29T19:41:00Z     2005-01-30T02:41:00Z     0               12548017        0        3.480471        Inf             
%csa_caa_irfu_test_core('WHISPER',test_location,3,1,'2005-01-29T19:41:00Z','2005-01-30T02:41:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('WHISPER',test_location,3,1,'2005-01-29T16:39:00Z','2005-01-29T20:39:00Z',fid);
file_size(4)=temp1;time_ratio(4)=temp2;

%C3_CP_STA_DWF_NBR                       2005-01-29T22:02:00Z     2005-01-30T08:02:00Z     123217894       123217809       22.393946       21.378309       0.954647        
%csa_caa_irfu_test_core('STAFF',test_location,3,1,'2005-01-29T22:02:00Z','2005-01-30T08:02:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('STAFF',test_location,3,1,'2005-01-29T03:59:00Z','2005-01-29T07:59:00Z',fid);
file_size(5)=temp1;time_ratio(5)=temp2;

%C2_CP_RAP_EPITCH                        2005-01-29T09:45:00Z     2005-01-30T04:45:00Z     19230808        19230861        24.375605       25.287603       1.037414        
%csa_caa_irfu_test_core('RAPID',test_location,2,2,'2005-01-29T09:45:00Z','2005-01-30T04:45:00Z',fid);
 [temp1,temp2]=csa_caa_irfu_test_core('RAPID',test_location,2,2,'2005-01-29T00:33:00Z','2005-01-29T21:33:00Z',fid);
file_size(6)=temp1;time_ratio(6)=temp2;

%C3_CP_RAP_DE                            2005-01-29T04:29:00Z     2005-01-29T14:29:00Z     5439168         5438771         11.185409       15.322274       1.369845        
%csa_caa_irfu_test_core('RAPID',test_location,3,3,'2005-01-29T04:29:00Z','2005-01-29T14:29:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('RAPID',test_location,3,3,'2005-01-29T16:11:00Z','2005-01-30T00:11:00Z',fid);
file_size(7)=temp1;time_ratio(7)=temp2;

%C3_CP_PEA_PADMARH_DEFlux                2005-01-29T15:42:00Z     2005-01-30T09:42:00Z     80358078        80357985        41.064095       102.186268      2.488458        
%csa_caa_irfu_test_core('PEACE',test_location,3,1,'2005-01-29T15:42:00Z','2005-01-30T09:42:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('PEACE',test_location,3,1,'2005-01-29T11:58:00Z','2005-01-29T14:58:00Z',fid);
file_size(8)=temp1;time_ratio(8)=temp2;

%C3_CP_FGM_FULL                          2005-01-29T06:40:00Z     2005-01-29T21:40:00Z     58121849        58121906        36.545009       56.645684       1.550025        
%csa_caa_irfu_test_core('FGM',test_location,3,1,'2005-01-29T06:40:00Z','2005-01-29T21:40:00Z',fid);
 [temp1,temp2]=csa_caa_irfu_test_core('FGM',test_location,3,1,'2005-01-29T20:38:00Z','2005-01-30T05:38:00Z',fid);
file_size(9)=temp1;time_ratio(9)=temp2;

%C3_CP_EFW_L2_V3D_GSE_EX                 2005-01-29T03:07:00Z     2005-01-29T14:07:00Z     59440347        0               44.729541       0        0        
%csa_caa_irfu_test_core('EFW',test_location,3,2,'2005-01-29T03:07:00Z','2005-01-29T14:07:00Z',fid);
%[temp1,temp2]=csa_caa_irfu_test_core('EFW',test_location,3,2,'2005-01-29T04:25:00Z','2005-01-29T15:25:00Z',fid);
%file_size(10)=temp1;time_ratio(10)=temp2;

%C1_CP_CIS-CODIF_HS_H1_PEF               2005-01-29T23:20:00Z     2005-01-30T13:20:00Z     0               0               0        0        NaN             
%csa_caa_irfu_test_core('CIS',test_location,1,1,'2005-01-29T23:20:00Z','2005-01-30T13:20:00Z',fid);

%C3_CP_CIS-CODIF_HS_H1_MOMENTS           2005-01-29T05:45:00Z     2005-01-29T11:45:00Z     140930          140861          7.092362        18.483133       2.606062        
%csa_caa_irfu_test_core('CIS',test_location,3,2,'2005-01-29T05:45:00Z','2005-01-29T11:45:00Z',fid);
 [temp1,temp2]=csa_caa_irfu_test_core('CIS',test_location,3,2,'2005-01-29T02:35:00Z','2005-01-29T07:35:00Z',fid);
file_size(11)=temp1;time_ratio(11)=temp2;

%C4_CP_CIS-CODIF_RPA_HE1_PSD             2005-01-29T12:41:00Z     2005-01-30T09:41:00Z     50486           50417           10.273389       3.749539        0.364976        
%csa_caa_irfu_test_core('CIS',test_location,4,3,'2005-01-29T09:41:00Z','2005-01-30T09:41:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('CIS',test_location,4,3,'2005-01-29T09:34:00Z','2005-01-29T15:34:00Z',fid);
file_size(12)=temp1;time_ratio(12)=temp2;

%C2_CP_PEA_PITCH_FULL_PSD                2005-01-29T23:32:00Z     2005-01-30T02:32:00Z     19596927        19596834        19.100827       36.705022       1.921646        
%csa_caa_irfu_test_core('PEACE',test_location,2,3,'2005-01-29T23:32:00Z','2005-01-30T02:32:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('PEACE',test_location,2,3,'2005-01-29T06:37:00Z','2005-01-29T12:37:00Z',fid);
file_size(13)=temp1;time_ratio(13)=temp2;

%C3_CP_WHI_ELECTRON_GYROFREQUENCY        2005-01-29T03:15:00Z     2005-01-29T23:15:00Z     33883           33792           10.175360       3.785975        0.372073        
%csa_caa_irfu_test_core('WHISPER',test_location,3,2,'2005-01-29T03:15:00Z','2005-01-29T23:15:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('WHISPER',test_location,3,2,'2005-01-29T19:58:00Z','2005-01-30T12:58:00Z',fid);
file_size(14)=temp1;time_ratio(14)=temp2;

%C3_CP_STA_PSD                           2005-01-29T06:48:00Z     2005-01-29T11:48:00Z     9914360         50829           14.357565       7.009980        0.488243        
%csa_caa_irfu_test_core('STAFF',test_location,3,2,'2005-01-29T06:48:00Z','2005-01-29T11:48:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('STAFF',test_location,3,2,'2005-01-29T08:35:00Z','2005-01-29T10:35:00Z',fid);
file_size(15)=temp1;time_ratio(15)=temp2;

%C3_CP_ASP_IONC                          2005-01-29T22:20:00Z     2005-01-30T02:20:00Z     261484          261407          11.017326       11.534812       1.046970        
%csa_caa_irfu_test_core('ASPOC',test_location,3,1,'2005-01-29T22:20:00Z','2005-01-30T02:20:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('ASPOC',test_location,3,1,'2005-01-29T21:52:00Z','2005-01-30T16:52:00Z',fid);
file_size(16)=temp1;time_ratio(16)=temp2;

%C1_CP_DWP_CORR_ST                       2005-01-29T06:36:00Z     2005-01-29T17:36:00Z     916009          916057          8.624390        10.161210       1.178195        
[temp1,temp2]=csa_caa_irfu_test_core('DWP',test_location,1,1,'2005-01-29T06:36:00Z','2005-01-29T17:36:00Z',fid);
file_size(17)=temp1;time_ratio(17)=temp2;

%C4_CG_WBD_GIFPLOT                       2005-01-29T08:49:00Z     2005-01-29T22:49:00Z     166190          166190          193.225015      295.171368      1.527604        
[temp1,temp2]=csa_caa_irfu_test_core('WBD',test_location,4,1,'2005-01-29T08:49:00Z','2005-01-29T22:49:00Z',fid);
file_size(18)=temp1;time_ratio(18)=temp2;

%C1_CP_EDI_EGD                           2005-01-29T13:55:00Z     2005-01-29T19:55:00Z     839810          839739          10.534330       50.271011       4.772113        
%csa_caa_irfu_test_core('EDI',test_location,1,1,'2005-01-29T13:55:00Z','2005-01-29T19:55:00Z',fid);
[temp1,temp2]=csa_caa_irfu_test_core('EDI',test_location,1,1,'2005-01-29T10:18:00Z','2005-01-29T13:18:00Z',fid);
file_size(19)=temp1;time_ratio(19)=temp2;

% %% Generate Figure and PNG of download time ratio CAA/CSA
figure
file_size(file_size==0)=NaN;
time_ratio(time_ratio==0)=NaN;
plot(file_size/1024^2,time_ratio,'*')

title(sprintf(['Download time ratio CAA/CSA \n',...
    ' IRFU test: ','\n'...
    'performed from ',num2str(test_location),' on ',num2str(date)]))
ylabel(sprintf(['t_{caa}/t_{csa}\n','CSA faster than CAA (x times)\n when ratio > 1']))
xlabel('File size [Mbytes]')
set(gca,'ylim',[0 min(9,max(time_ratio))+1]);
xlims=get(gca,'xlim'); xlims(1)=0; set(gca,'xlim',xlims); line(xlims,[1 1],'color','k');
saveas(gcf,['IRF_test_from_',num2str(test_location),'_',num2str(date)],'png')

fclose(fid);
