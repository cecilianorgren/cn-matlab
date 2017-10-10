% script_download_ib

cd /Users/Cecilia/Data/IB;

t1=textread('/Users/Cecilia/Data/IB/IB_list_JK.txt','%s %*[^\n]');
l=0;p=0;
for k=1:length(t1)
    switch t1{k}(17)
        case '3'
            l=l+1;
            year(l,1)=str2num(['20', t1{k}(1:2)]);
            month(l,1)=str2num(t1{k}(3:4));
            day(l,1)=str2num(t1{k}(5:6));
            hour(l,1)=str2num(t1{k}(7:8));
            minute(l,1)=str2num(t1{k}(9:10));
            second(l,1)=str2num(t1{k}(11:12));
            str3{l,:}=['20',t1{k}(1:6)];
        case '4'
            p=p+1;
            year(p,2)=str2num(['20', t1{k}(1:2)]);
            month(p,2)=str2num(t1{k}(3:4));
            day(p,2)=str2num(t1{k}(5:6));
            hour(p,2)=str2num(t1{k}(7:8));
            minute(p,2)=str2num(t1{k}(9:10));
            second(p,2)=str2num(t1{k}(11:12));
            str4{p,:}=['20',t1{k}(1:6)];
    end
end
tint3=[toepoch([year(:,1) month(:,1) day(:,1) hour(:,1) minute(:,1) second(:,1)])-2,...
       toepoch([year(:,1) month(:,1) day(:,1) hour(:,1) minute(:,1) second(:,1)])+14];
tint4=[toepoch([year(:,2) month(:,2) day(:,2) hour(:,2) minute(:,2) second(:,2)])-2,...
       toepoch([year(:,2) month(:,2) day(:,2) hour(:,2) minute(:,2) second(:,2)])+14];
clear year month day hour minute second k l p t1
%%
for k=1:size(tint3,1)
    str=[num2str(year(k,1)),num2str(month(k,1)),num2str(day(k,1))];
    eval(['mkdir ',str3])
end
%%
for k=1:size(tint3,1)
    %str=[num2str(year(k,1)),num2str(month(k,1)),num2str(day(k,1))];
    eval(['cd /Users/Cecilia/Data/IB/',str3{k}]);   
    caa_download(tint3(k,:),'C3_CP_EFW_L2_EB')
    caa_download(tint4(k,:),'C4_CP_EFW_L2_EB')    
end

%% Download Magnetic field data for a longer period
tint_long=[tint3(:,1)-600 tint3(:,2)+600]; % +- 10 min

for k=2:size(tint_long,1)    
    eval(['cd /Users/Cecilia/Data/IB/',str3{k}]);   
    caa_download(tint_long(k,:),'C3_CP_FGM_FULL_ISR2')
    caa_download(tint_long(k,:),'C4_CP_FGM_FULL_ISR2')
    caa_download(tint_long(k,:),'C3_CP_EFW_L2_E3D_INERT')
    caa_download(tint_long(k,:),'C4_CP_EFW_L2_E3D_INERT')
end
