function [file_size,time_ratio,time_range]=csa_caa_irfu_test_core(experiment,test_location,satnumber,filename_number,starts,ends,fid)
csa_server='http://csa.esac.esa.int/csa/aio/product-action?';
uname_csa='USERNAME=amasson&';
pass_csa='PASSWORD=Asdfgh(9&';
caa_server='http://caa.estec.esa.int/caa_query?';
uname_caa='uname=Masson&';
pass_caa='pwd=ob900i&';

switch experiment
    case 'ASPOC'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  ASPOC %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Science
        %%
        filenames{1}='CP_ASP_IONC';
        
        
        
    case 'CIS'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  CIS   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Science
        %% Moments
        filenames{1}='CP_CIS-CODIF_HS_H1_PEF';
        filenames{2}='CP_CIS-CODIF_HS_H1_MOMENTS';
        filenames{3}='CP_CIS-CODIF_RPA_He1_PSD';
        
        
    case 'DWP'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  DWP   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Science
        %%
        filenames{1}='CP_DWP_CORR_ST';
        
    case 'EDI'
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  EDI   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        filenames{1}='CP_EDI_EGD';
        
    case 'EFW'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  EFW   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        filenames{1}='CP_EFW_L3_E3D_INERT';
        filenames{2}='CP_EFW_L2_V3D_GSE_EX';
        
    case 'FGM'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  FGM   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        filenames{1}='CP_FGM_FULL';
        
    case 'PEACE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  PEACE   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        filenames{1}='CP_PEA_3DXPH_cnts';
        filenames{2}='CP_PEA_PADMARH_DEFlux';
        filenames{3}='CP_PEA_PITCH_FULL_PSD';
        
    case 'RAPID'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  RAPID   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        filenames{1}='CP_RAP_PAD_L3DD';
        filenames{2}='CP_RAP_EPITCH';
        filenames{3}='CP_RAP_DE';
        
        
    case 'STAFF'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  STAFF   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filenames{1}='CP_STA_DWF_NBR';
        filenames{2}='CP_STA_PSD';
        
    case 'WBD'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  WBD   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filenames{1}='CG_WBD_GIFPLOT';
        
        
    case 'WHISPER'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  WHISPER   test %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Science
        %%
        filenames{1}='CP_WHI_ELECTRON_DENSITY';
        filenames{2}='CP_WHI_ELECTRON_GYROFREQUENCY';
        
end

time_range=['time_range=',num2str(starts),'/',num2str(ends),'&refdoc=0'];
% log_filename=['IRF-U_test_from_',num2str(test_location),'_',num2str(date),'.txt'];
% fid=fopen(num2str(log_filename),'w');

%% Download specific dataset
for j=filename_number:filename_number
    
    for satnum=satnumber:satnumber
        
        filename=['dataset_id=C',num2str(satnum),'_',filenames{j},'&'];
        URL_CAA=[caa_server,uname_caa,pass_caa,filename,time_range,'&nonotify=1'];
        CAAfileName=tempname;
        gz_CAAFileName_temp = [CAAfileName '.zip'];tic
        gz_CAAFileName=urlwrite(URL_CAA,gz_CAAFileName_temp);
        caa_download_time=toc;
        [~,test_empty]=grep('-s',{'returned'},gz_CAAFileName);
        
        if isempty(strcmp(test_empty.match,'***CAA-Error*** Your request returned no results'))
            %           fileNames=unzip(gzFileName);
            
            if ~isempty(gz_CAAFileName)
                file_caa=dir(gz_CAAFileName);
            else
                file_caa.bytes=0;
            end
        else
            file_caa.bytes=0;
        end
        
        URL_CSA=[csa_server,uname_csa,pass_csa,filename,'START_DATE=',starts,'&END_DATE=',ends,'&NON_BROWSER'];
        CSAfileName=tempname;
        gz_CSAFileName_temp = [CSAfileName '.gz'];tic
        [gz_CSAFileName,csa_status]=urlwrite(URL_CSA,gz_CSAFileName_temp);
        csa_download_time=toc;
        
        if csa_status==1
            
            file_csa=dir(gz_CSAFileName);
        else
            file_csa.bytes=0;
        end
        
        %% Display the end result of the comparative download in the command window and log it
        if  file_csa.bytes==0 && file_caa.bytes ==0
            %             cprintf('comment',[num2str(filename(12:end-1)),'\tCSA: ',num2str(file_csa.bytes),' bytes, ',...
            %                 num2str(csa_download_time),'s; CAA: ',num2str(file_caa.bytes),' bytes, ',...
            %                 num2str(caa_download_time),'s\n']);
            %   elseif  abs(file_caa.bytes-file_csa.bytes)>(file_caa.bytes/100)&&abs(file_caa.bytes-file_csa.bytes)>2000
            %% Test precise content
        elseif  abs(file_caa.bytes-file_csa.bytes)>(file_caa.bytes/100)&&abs(file_caa.bytes-file_csa.bytes)>2000
            %%
            test_outcome_message=[num2str(filename(12:end-1)),'       CSA: ',num2str(file_csa.bytes),' bytes, ',...
                num2str(csa_download_time),'s; CAA: ',num2str(file_caa.bytes),' bytes, ',...
                num2str(caa_download_time),'s\n'];
            cprintf('-err',test_outcome_message);
            if file_caa.bytes > 0 && file_csa.bytes > 0
                time_ratio(j+satnum*length(filenames))=caa_download_time/csa_download_time;
                file_size(j+satnum*length(filenames))=file_csa.bytes;
            end
            fprintf(fid,'%s \r\n',test_outcome_message(1:end-2));
            if file_caa.bytes > 0 && file_csa.bytes > 0
                CAAfilenames=unzip(gz_CAAFileName); NUMBER_of_FILES=size(CAAfilenames,2);CAA_downloadsize=0;
                for i=1:NUMBER_of_FILES,temp=dir(CAAfilenames{i});CAA_downloadsize=temp.bytes+CAA_downloadsize;end
                gunzip(gz_CSAFileName);
                CSAfilenames=untar(CSAfileName); NUMBER_of_FILES=size(CSAfilenames,2);CSA_downloadsize=0;
                for i=1:NUMBER_of_FILES,temp2=dir(CSAfilenames{i});CSA_downloadsize=temp2.bytes+CSA_downloadsize;end
                diffsize_perc_unzip =  abs(CSA_downloadsize-CAA_downloadsize)/CAA_downloadsize;
                test_outcome_message2=[num2str(filename(12:end-1)),'      Uncompressed CSA: ',num2str(CSA_downloadsize),' bytes, ',...
                    'CAA: ',num2str(CAA_downloadsize),' bytes, Difference: ',num2str(diffsize_perc_unzip*100),' percent\n'];
                if diffsize_perc_unzip<0.1, cprintf('comment',test_outcome_message2); else cprintf('-err',test_outcome_message2); end;
                
                fprintf(fid,'%s \r\n',test_outcome_message2(1:end-2));
            end
        elseif caa_download_time<(csa_download_time/2)
            test_outcome_message= [num2str(filename(12:min(34,end-1))),'\t\t',num2str(starts),'\t',num2str(ends),'\t',num2str(file_caa.bytes),'\t',...
                num2str(file_csa.bytes),' \t',num2str(caa_download_time),' \t',...
                num2str(csa_download_time),'\t',...
                num2str(csa_download_time/caa_download_time),'\n'];
            cprintf([1 0.5 0], test_outcome_message);
            time_ratio(satnum+4*(j-1))=caa_download_time/csa_download_time;
            file_size(satnum+4*(j-1))=file_csa.bytes;
            fprintf(fid,'%s \r\n',test_outcome_message(1:end-2));
        else
            CAAfilenames=unzip(gz_CAAFileName); NUMBER_of_FILES=size(CAAfilenames,2);CAA_downloadsize=0;
            for i=1:NUMBER_of_FILES,temp=dir(CAAfilenames{i});CAA_downloadsize=temp.bytes+CAA_downloadsize;end
            gunzip(gz_CSAFileName);
            CSAfilenames=untar(CSAfileName); NUMBER_of_FILES=size(CSAfilenames,2);CSA_downloadsize=0;
            for i=1:NUMBER_of_FILES,temp2=dir(CSAfilenames{i});CSA_downloadsize=temp2.bytes+CSA_downloadsize;end
            diffsize_perc_unzip = 100* abs(CSA_downloadsize-CAA_downloadsize)/CAA_downloadsize;
            %              test_outcome_message= [num2str(filename(12:min(34,end-1))),'\t\t',num2str(starts),'\t',num2str(ends),'\t',num2str(file_caa.bytes),'\t',...
            %                 num2str(file_csa.bytes),' \t',num2str(caa_download_time),' \t',...
            %                 num2str(csa_download_time),'\t',...
            %                 num2str(csa_download_time/caa_download_time),'\n'];
            if strcmp(filename(12:min(34,end-1)),'C3_CP_RAP_DE')
                test_outcome_message= [num2str(filename(12:min(34,end-1))),'\t\t\t',num2str(starts),'\t',num2str(ends),'\t',num2str(CAA_downloadsize),'\t',...
                    num2str(CSA_downloadsize),' \t',num2str(caa_download_time),' \t',...
                    num2str(csa_download_time),'\t',...
                    num2str(csa_download_time/caa_download_time),'\t',num2str(diffsize_perc_unzip),'\n'];
                 fprintf(fid,'%s\t\t\t %s\t %s\t %d\t %d\t %f\t %f\t %f\t %f\t \r\n',filename(12:min(34,end-1)),starts,ends,CAA_downloadsize,...
                    CSA_downloadsize,caa_download_time,csa_download_time,csa_download_time/caa_download_time,diffsize_perc_unzip);
            else
                test_outcome_message= [num2str(filename(12:min(34,end-1))),'\t\t',num2str(starts),'\t',num2str(ends),'\t',num2str(CAA_downloadsize),'\t',...
                    num2str(CSA_downloadsize),' \t',num2str(caa_download_time),' \t',...
                    num2str(csa_download_time),'\t',...
                    num2str(csa_download_time/caa_download_time),'\t',num2str(diffsize_perc_unzip),'\n'];
                    fprintf(fid,'%s\t\t %s\t %s\t %d\t %d\t %f\t %f\t %f\t %f\t \r\n',filename(12:min(34,end-1)),starts,ends,CAA_downloadsize,...
                    CSA_downloadsize,caa_download_time,csa_download_time,csa_download_time/caa_download_time,diffsize_perc_unzip);
            end
            cprintf(test_outcome_message);
            if file_caa.bytes > 0 && file_csa.bytes > 0
                time_ratio(satnum+4*(j-1))=caa_download_time/csa_download_time;
                file_size(satnum+4*(j-1))=file_csa.bytes;
            end
            %fprintf(fid,'%s \n',sprintf(test_outcome_message));
        end
        %         if file_csa.bytes>0&&file_caa.bytes>0
        %             if ~strcmp(file_caa.name,file_csa(i).name)
        %                 cprintf(['Downloaded fileame is different from CSA ',num2str(file_csa(i).name),...
        %                     ' wrt CAA\n ',num2str(file_caa.name),'\n']);
        %             end
        %         end
        
    end
end
if ~isempty(CSA_downloadsize)
    file_size=CSA_downloadsize;
else
    file_size=0;
end

time_ratio=caa_download_time/csa_download_time;

