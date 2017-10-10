gsm0E3=c_coord_trans('dsi','gsm',[diE3(:,1:3) diE3(:,1)*0],'cl_id',3);
gsm0E4=c_coord_trans('dsi','gsm',[diE4(:,1:3) diE4(:,1)*0],'cl_id',4);

gsmzE3=c_coord_trans('dsi','gsm',diE3,'cl_id',3);
gsmzE4=c_coord_trans('dsi','gsm',diE4,'cl_id',4);
%%
facgsm0E3=cn_m_trans(cn_toepoch(t1,t2,gsm0E3),gsmM,1);
facgsm0E4=cn_m_trans(cn_toepoch(t1,t2,gsm0E4),gsmM,1);

facgsmzE3=cn_m_trans(cn_toepoch(t1,t2,gsmzE3),gsmM,1);
facgsmzE4=cn_m_trans(cn_toepoch(t1,t2,gsmzE4),gsmM,1);
%%
figure;irf_plot({facgsm0E3,facgsmzE3},'comp')
userdata=get(gcf,'userdata');
legend('diEz=0','EdotB=0')
title(userdata.subplot_handles(1),'E FAC C3')
set(gcf,'PaperPositionMode','auto');
figure;irf_plot({facgsm0E4,facgsmzE4},'comp')
userdata=get(gcf,'userdata');
legend('diEz=0','EdotB=0')
title(userdata.subplot_handles(1),'E FAC C4')
set(gcf,'PaperPositionMode','auto');
%%
figure;irf_plot({cn_toepoch(t1,t2,gsm0E3),cn_toepoch(t1,t2,facgsmzE3)},'comp')
userdata=get(gcf,'userdata');
legend('diEz=0','EdotB=0')
title(userdata.subplot_handles(1),'E FAC C3')
set(gcf,'PaperPositionMode','auto');
figure;irf_plot({facgsm0E4,facgsmzE4},'comp')
userdata=get(gcf,'userdata');
legend('diEz=0','EdotB=0')
title(userdata.subplot_handles(1),'E FAC C4')
set(gcf,'PaperPositionMode','auto');

%% B parallel
