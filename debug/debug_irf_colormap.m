close all; figure(71); h(1)=subplot(3,1,1); h(2)=subplot(3,1,2); h(3)=subplot(3,1,3);%axes;
axis_handle1 = irf_match_phibe_vis('velocity/density',h(1),n,v,corr_v);
axis_handle2 = irf_match_phibe_vis('velocity/density',h(2),n,v,corr_v);
axis_handle3 = irf_match_phibe_vis('velocity/density',h(3),n,v,corr_v);
cStrPrint = [irf_time(tint(1),'yyyymmddhhmmss') '-' irf_time(tint(2),'yyyymmddhhmmss') '_f' num2str(f_highpass) 'Hz_sc' num2str(sc) '_c'];
%print('-dpng',['corr3' '.png'])
irf_colormap(axis_handle2.ax,'default')

    
%%
allall = findall(gcf);
texts = []; 
for jj=1:numel(allall); 
    type = get(allall(jj),'type'); 
    %disp([get(allall(jj),'type')]); 
    if 1
    if strcmp(type,'text') 
        alltexts = [texts allall(jj)]; 
        disp([num2str(jj) ' ' get(allall(jj),'string')]); 
    end 
    end
end
%%
allobj = findobj(gcf);
objtexts = [];
for jj=1:numel(allobj); 
    type = get(allobj(jj),'type'); 
    disp([get(allobj(jj),'type')]); 
    if strcmp(type,'text') 
        objtexts = [texts allobj(jj)]; 
        disp([num2str(jj) ' ' get(allobj(jj),'string')]); 
    end 
end

%%
allaxes = findall(gcf,'type','axes');
allchildren = get(allaxes,'Children');

for gg = 1:numel(allchildren)
    get(allchildren{gg})
end

%%
close all; figure(71); 
for kk=1:3; h(kk)=subplot(3,1,kk); peaks; colorbar; end

cn.colormap(h(2),'space')






