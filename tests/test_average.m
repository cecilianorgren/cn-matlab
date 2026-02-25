if 0
for k=1:size(gsmE3,1)/900
    k
        gsmE3lowres(k,1)=mean(gsmE3((1:900)*k),1);
        gsmE3lowres(k,2)=mean(gsmE3((1:900)*k),2);
end
end
%% E with diEz=0;
E3=c_coord_trans('dsi','gsm',[diE3(:,1:3) diE3(:,1)*0],'cl_id',3);
%%
original=E3;
freq=450;
dim=size(original,1);
rest=mod(dim,freq);



dim_1=tocolumn(nanmean(reshape(original(1:(dim-rest),1),freq,(dim-rest)/freq),1));
dim_2=tocolumn(nanmean(reshape(original(1:(dim-rest),2),freq,(dim-rest)/freq),1));
dim_3=tocolumn(nanmean(reshape(original(1:(dim-rest),3),freq,(dim-rest)/freq),1));
dim_4=tocolumn(nanmean(reshape(original(1:(dim-rest),4),freq,(dim-rest)/freq),1));
E3lowres=[dim_1 dim_2 dim_3 dim_4];
%%
figure;irf_plot(E3(:,[1 3]),'r');hold(gca,'on');irf_plot(E3lowres(:,[1 3]),'k');
figure;irf_plot(gsmE3(:,[1 3]),'r');hold(gca,'on');irf_plot(gsmE3lowres(:,[1 3]),'k');
figure;irf_plot({gsmE3, E3},'comp');
legend(gca,'EdotB=0','diEz=0')

%%
absdiEz0=irf_abs(E3);
absEdotB=irf_abs(gsmE3);
figure;irf_plot({absEdotB,absdiEz0},'comp')
legend(gca,'EdotB=0','diEz=0')
%%
figure;irf_plot({c_coord_trans('dsi','gsm',diE3(:,1:4),'cl_id',3),...
    c_coord_trans('dsi','gsm',[diE3(:,1:3) diE3(:,1)*0],'cl_id',3)},'comp')
legend(gca,'EdotB=0','diEz=0')