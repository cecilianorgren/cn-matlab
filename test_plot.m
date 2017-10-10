h=cn_plot3d([0 0 0],[0 0 0],' ');
axis equal; set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
for k=1:100:542
    cn_plot3d(h,[0 0 0],b3norm(k,:),' ',b4norm(k,:),' ')
end

%%
v=500; % km/s
freq=450;
dz=v/freq;
step=5;
scale=10;
h=cn_plot3d([-5 -5 30],[0 0 0],'b3');
cn_plot3d(h,[5 5 0],[0 0 0],'b4');
par=0;
for k=1:step:size(b3,1)/3   
    cn_plot3d(h,[-5 -5 par+30],facE3(k,2:4)./([1 1 1]*scale),'');
    cn_plot3d(h,[5 5 par],facE4(k,2:4)./([1 1 1]*scale),'');
    par=par-dz*step;
end
axis equal

%%
kdir=kdir3;
corrx=corrx3(:,1);
corrmax=max(corrx);
maxi=find(corrx==corrmax);
corrmin=min(corrx);
%%
corrmax=1;
s=1;
figure;
zhat=irf_norm(mean(cn_toepoch(t1,t2,gsmB3),1));
h=cn_plot3d([0 0 0],-x*corrmax*s,'MVAx',propx,'old tool',zhat(2:4),'B');
str='';
for k=1:size(kdir,1)
    if k==maxi;
        str='max';
    end
    cn_plot3d(h,[0 0 0],s*corrx(k)*kdir(k,:),str);
    str='';
end
axis equal; set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
%%
for k=1:size(xc,1)
    cn_plot3d(h,[0 0 0],xc(k,:),str);
end
