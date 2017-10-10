if ~exist('gsmE3','var'); load matlabE; end
if ~exist('gsmB3','var'); load matlabB; end
% t1=[2007 08 31 10 19 06.90]; t2=[2007 08 31 10 19 07.50]; % lhdw
t1=[2007 08 31 10 18 42.45]; t2=[2007 08 31 10 18 43.55]; % turb + dl + space 
t1=[2007 08 31 10 18 42.50]; t2=[2007 08 31 10 18 42.85]; % turb
t1=[2007 08 31 10 18 42.50]; t2=[2007 08 31 10 18 43.00]; % turb + dl
t1=[2007 08 31 10 18 42.70]; t2=[2007 08 31 10 18 43.00]; % dl - most turb
tint=[toepoch(t1) toepoch(t2)];
sclist=4;
%% Make with smaller E, only the measured one
%c_eval('size([diE?(:,1:3) zeros(size(diE?,1),1)])',sclist)
c_eval('diE?noz=irf_tlim(diE?(:,1:3),[tint(1)-10 tint(2)+10]);',sclist)
c_eval('diE?noz=[diE?noz(:,1:3) diE3noz(:,1)*0];',sclist)
c_eval('gsmE?noz=c_coord_trans(''dsi'',''gsm'',diE?noz,''cl_id'',?);',sclist)
%% make parallel
c_eval('diE?par=irf_edb(diE?noz,diB?,70,''Epar'');',sclist);
c_eval('gsmE?par=c_coord_trans(''dsi'',''gsm'',diE?par,''cl_id'',?);',sclist)
%% Choose which spacecraft to use
sc=4;
c_eval('gsmB=gsmB?fgm;',sc);
c_eval('gsmE=gsmE?;',sc);

% Direction
f_filt=1; % Hz
% with cn_toepoch [x y z corr phiE Bz Ek En ufEn ufEk]=tool_direction(cn_toepoch(t1,t2,gsmB),cn_toepoch(t1,t2,gsmE),300,f_filt);
% with irf_tlim
[x y z corr phiE Bz Ek En ufEn ufEk]=tool_direction(irf_tlim(gsmB,tint),irf_tlim(gsmE,tint),300,f_filt);
index_dir=find(corr(:,1)==max(corr(:,1)));
x(index_dir,:)

%% Velocity
if size(gsmB,2)==4; gsmB=irf_abs(gsmB); end
B0=mean(irf_tlim(gsmB,[toepoch(t1) toepoch(t2)]),1); B0=B0(5);
vint=[100,2000];
title_str=['C',num2str(sc,'%0.f'),', GSM,   f_{filt} = ',num2str(f_filt,'%0.f'),...
    'Hz,  khat = [',num2str(x(index_dir,1),'%0.2f'),' ',...
    num2str(x(index_dir,2),'%0.2f'),' ',num2str(x(index_dir,3),'%0.2f'),']'];
n=linspace(0.01,0.1,10);
v=linspace(vint(1),vint(2),10);
[correl,vis]=tool_velocity(B0,Bz,phiE(:,[1 index_dir+1]),n,v,title_str,1);
index_v=find(correl(:,1)==min(correl(:,1)));
v(index_v)
tint=[Bz(1,1) Bz(end,1)];
%imwrite(im,map,['/Users/Cecilia/DF/Pics/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_C',num2str(sc,'%0.f'),'_f',num2str(sc,'%0.f'),'_velocity.gif'],'DelayTime',0.01,'LoopCount',inf);
%% trying to get 2d matrix
vint=[0.1,2490]; % km/s
nint=[0.05,1.0]; % cc
nv=600;
nn=550;
vs=linspace(vint(1),vint(2),nv);
%ns=linspace(nint(1),nint(2),nn);
%vs=logspace(log10(vint(1)),log10(vint(2)),nv);
ns=logspace(log10(nint(1)),log10(nint(2)),nn);
[correlation]=tool_velocity2(B0,Bz,phiE(:,[1 index_dir+1]),ns,vs,title_str);
figure(92);
pcolor(v,n,log10(correl));
p=gca;
hc=colorbar;


title(['Tool correlation matrix, amplitude,  khat = [',...
    num2str(x(index_dir,1),'%0.2f'),' ',...
    num2str(x(index_dir,2),'%0.2f'),' ',...
    num2str(x(index_dir,3),'%0.2f'),'] GSM'])
xlabel('Velocity [km/s]')
ylabel('Density [cc]')
ylabel(hc,'log10(sum((LHS-RHS)^2))')
[rows1,cols1,vals1]=find(correlation==repmat(min(correlation'),nv,1)');
[rows2,cols2,vals2]=find(correlation==repmat(min(correlation),nn,1));
hold on;
shading flat
plot(vs(cols1(1:(end-100))),ns(rows1(1:(end-100))),'color',[0 0 0],'linewidth',1)
plot(vs(cols2(1:(end-100))),ns(rows2(1:(end-100))),'color',[0 0 0],'linewidth',1)
set(gca,'yscale','log','TickDir','out')
hold off
%% Illustrate
title_str=['C',num2str(sc,'%0.f'),', GSM,   f_{filt} = ',num2str(f_filt,'%0.f'),'Hz'];
[A,im,map]=tool_gif(x,y,z,corr,[phiE(:,1) phiE(:,2:end)*v(index_v)],ufEk,ufEn,Bz,title_str);
tint=[Bz(1,1) Bz(end,1)];
imwrite(im,map,['/Users/Cecilia/DF/Pics/tool_T',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(2)),'HHMMSSFFF'),'_C',num2str(sc,'%0.f'),'_f',num2str(f_filt,'%0.f'),'_.gif'],'DelayTime',0.01,'LoopCount',inf);

%% Illustrate
if 0 
h=irf_plot(11);
for k=1:11
    index=index0+(k-6);
    irf_plot(h(k),{[phiE(:,1) phiE(:,index+1)./repmat(max(phiE(:,index+1)),size(phiE,1),1)],...
        [Bz(:,1) Bz(:,2)./repmat(max(Bz(:,2)),size(Bz,1),1)]},'comp')
    ylabel(h(k),[' ',num2str(index),' '])
end
end