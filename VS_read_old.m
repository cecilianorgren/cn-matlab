VirtualSatellite=textread('/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces1000.txt');
%%
N = 2 ; % every 2nd row
Myformat = ['%d%d' repmat('%*s',1,n-1)] ; 
% i.e., read the two columns as numbers and skip N-1 lines
[Col1, Col2] = textread('MyTextFile.txt', MyFormat, 'delimiter','\n')
%%
for k=1:27
    eval(['Point',num2str(k)])=VirtualSatellite;
end
%[t1,t2,tint_comments]=textread('/Users/Cecilia/Exjobb/BM/BM.txt','%s%s%[^\n]');
%% Clear variables
eval(['clear B',num2str(p),';'])
%% TEXTSCAN Bad example
fid=fopen('/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces1000.txt');
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
for p=1:27
    eval(['B',num2str(p),'=[];'])
    eval(['P',num2str(p),'=[Pos{1,1}(',num2str(p),',:) Pos{1,2}(',num2str(p),',:) Pos{1,3}(',num2str(p),',:)];']);
end
for k=1:10 % ten first time steps
    Data=textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
    for p=1:27 
        eval(['B',num2str(p),'=[B',num2str(p),'; Data{1,1}(',num2str(p),'), Data{1,2}(',num2str(p),'), Data{1,3}(',num2str(p),')];']);
    end
end
% Vector with the 27 positions, one xyz for each row
xyz=[Pos{1} Pos{2} Pos{3}];
%% TEXTSCAN Better example
format = repmat('%f',1,14*27); % 14*27 numbers
for core=1000
eval(['fid=fopen(''/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces',num2str(core),'.txt'');'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
Position=[Pos{1} Pos{2} Pos{3}];
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:378,length(Data)/27,27);
for k=1:27
    % Divides into 27 satellites
    eval(['B',num2str(k),' = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E',num2str(k),' = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Ji',num2str(k),' = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Je',num2str(k),' = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ni',num2str(k),' = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ne',num2str(k),' = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end
%% Plot some of the data to test
subplot(1,5,1)
plot(B1)
subplot(1,5,2)
plot(B2)
subplot(1,5,3)
plot(B3)
subplot(1,5,4)
plot(B4)
subplot(1,5,5)
plot(B5)
%% TEXTSCAN Yet better example
format = repmat('%f',1,14*27); % 14*27 numbers
core=[0 100 101 102 103 104 105 106 107 1000 1001 1002 1003 1004 1057 1058 1059 1060,...
    1061 1062 1063 1064 1065 1066 1074 1075 1076 1077 1078 1079 1080];
for p=1:2
eval(['fid=fopen(''/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces',num2str(core(p)),'.txt'');'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
eval(['Position',num2str(core(p)),'=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:378,length(Data)/27,27);
for k=1:27
    % Structure: QuantityCore(time,xyz,spot within core)
    eval(['B',num2str(core(p)),'(:,:,',num2str(k),') = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E',num2str(core(p)),'(:,:,',num2str(k),') = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Ji',num2str(core(p)),'(:,:,',num2str(k),') = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Je',num2str(core(p)),'(:,:,',num2str(k),') = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ni',num2str(core(p)),'(:,',num2str(k),') = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ne',num2str(core(p)),'(:,',num2str(k),') = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end
%% TEXTSCAN Even yet better example
format = repmat('%f',1,14*27); % 14*27 numbers
core=[0 100 101 102 103 104 105 106 107 1000 1001 1002 1003 1004 1057 1058 1059 1060,...
    1061 1062 1063 1064 1065 1066 1074 1075 1076 1077 1078 1079 1080];
for p=1:length(core)
eval(['fid=fopen(''/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces',num2str(core(p)),'.txt'');'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
eval(['Position(',num2str(p),',:,:)=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:length(Data),length(Data)/27,27);
for k=1:27
    % Structure: Quantity(core,sc (1-27),time,xyz)
    eval(['B(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Ji(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Je(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ni(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ne(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end
%% Plot some of the data to test
nplot=5;
for k=1:nplot
eval(['subplot(1,',num2str(nplot),')'])
subplot(1,5,1)
end
%% Searching of spots
xrange=[9 16]; xrange=[24 31];
yrange1=[5 6]; yrange2=[9 10];
zrange=[0 10];
%[c,sc]=ind2sub(size(Position),find(Position(:,:,1)>9));
[c1,sc1]=find(Position(:,:,1)>9);
[c2,sc2]=find(Position(:,:,1)<16);
[cA scA]=intersect([c1,sc1],[c2,sc2],'rows');
[c3,sc3]=find(Position(:,:,1)>24);
[c4,sc4]=find(Position(:,:,1)<31);
[cB scB]=intersect([c3,sc3],[c4,sc4],'rows');
[c5,sc5]=find(Position(:,:,2)>5);
[c6,sc6]=find(Position(:,:,2)>6);
[cC scC]=intersect([c5,sc5],[c6,sc6],'rows');
[c7,sc7]=find(Position(:,:,2)>9);
[c8,sc8]=find(Position(:,:,2)<10);
[cD scD]=intersect([c7,sc7],[c8,sc8],'rows');
%% Searching of spots
xrange=[9 16]; xrange=[24 31];
yrange1=[5 6]; yrange2=[9 10];
zrange=[0 10];
%[c,sc]=ind2sub(size(Position),find(Position(:,:,1)>9));
[c1,sc1]=find(Position(:,:,1)>11);
[c2,sc2]=find(Position(:,:,1)<15);
csc1=intersect([c1,sc1],[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,2)>8);
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,2)<11);
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,3)>4);
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,3)<7);
csc1=intersect(csc1,[c2,sc2],'rows');

csc1

%% Plot some of the data that was found
m=1;
plot(E(csc1(m,1),csc1(m,2),:,1))
%%
[c2,sc2]=find(Position(:,:,1)<31);
[c1 sc1]=intersect([c1,sc1],[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,1)>24);
[c1 sc1]=intersect([c1,sc1],[c2,sc2],'rows');

[c2,sc2]=find(Position(:,:,2)>9);
[c1 sc1]=intersect([c1,sc1],[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,2)<10);
[c1 sc1]=intersect([c1,sc1],[c2,sc2],'rows');
%% Check distributions of satellites
binsize=1;
matrix=zeros(40,15);
for x=0:39
    for y=0:14
        hej=tred_locate([x x+1],[y y+1],[0 10],Position);
        matrix(x+1,y+1)=size(hej,1);
    end
end
surf(repmat(0:40,16,1)',repmat(0:15,41,1),zeros(41,16),matrix); 
view([0 0 1]); colorbar; xlabel('x'); ylabel('y'); title('Number of spots summed over z')
%% Plot fields
figure;

plot(squeeze(E(find(core==1001),1,:,1)),'b');hold on;
plot(squeeze(E(find(core==1001),1,:,2)),'g')
plot(squeeze(E(find(core==1001),1,:,3)),'r')
%%
figure;
x=1:size(E,3);
if 0 % Closelying cores
plot(x,squeeze(E(find(core==1000),1,:,2)),...
     x,squeeze(E(find(core==1001),1,:,2)),...
     x,squeeze(E(find(core==1002),1,:,2)),...
     x,squeeze(E(find(core==1003),1,:,2)),...
     x,squeeze(E(find(core==1004),1,:,2)))
end
subplot(3,1,1)
if 1 % 1 core different satellites
plot(x,squeeze(B(find(core==1000),1,:,1)),...
     x,squeeze(B(find(core==1000),2,:,1)),...
     x,squeeze(B(find(core==1000),3,:,1)),...
     x,squeeze(B(find(core==1000),4,:,1)),...
     x,squeeze(B(find(core==1000),5,:,1)))
end
subplot(3,1,2)
if 1 % 1 core different satellites
plot(x,squeeze(E(find(core==1000),1,:,2)),...
     x,squeeze(E(find(core==1000),2,:,2)),...
     x,squeeze(E(find(core==1000),3,:,2)),...
     x,squeeze(E(find(core==1000),4,:,2)),...
     x,squeeze(E(find(core==1000),5,:,2)))
end
subplot(3,1,3)
plot(x,squeeze(E(find(core==1000),4,:,2))*1e0,...
     x,squeeze(B(find(core==1000),4,:,2)))
 
%% Reading all electric fields
x=(1:length(E(1,1,:,1)));
for k=1:30
    for p=1:9:27
        figure;
        set(gcf,'PaperPositionMode','auto',...
            'Position',[1     1   471   927]);
        h=irf_plot(9);
        for l=1:9
            irf_plot(h(l),[tocolumn(x) squeeze(squeeze(Epar(k,p+l-1,:,:)))]);
            %subplot(9,1,l)
            %plot(x,squeeze(E(k,p+l-1,:,1)),x,squeeze(E(k,p+l-1,:,2)))
            if l~=9; set(gca,'XTick',[]); end                
        end
        eval(['print -dpng /Users/Cecilia/VTK/Satellitbilder/Epar_c',num2str(core(k)),'sc',num2str(p),'-',num2str(p+8),'.png']);
    end
    close all;
end
%% Calculate Epar
x=(1:length(E(1,1,:,1)));
for k=1:30
    for l=1:27
        Epar(k,l,:,1)=irf_dot(squeeze(E(k,l,:,:)),irf_norm(squeeze(B(k,l,:,:))));
    end
end
%% Plot position of cores
figure;
for k=1:30
    for l=14
        plot3(Position(k,l,1),Position(k,l,2),Position(k,l,3)); hold on
        text(Position(k,l,1),Position(k,l,2),Position(k,l,3),num2str(core(k)))
    end
end
set(gca,'XLim',[0 40],'YLim',[0 15],'ZLim',[0 10])
xlabel('x')
ylabel('y')
zlabel('z')

%% Find position
fid=fopen('/Users/Cecilia/TRED46/DATA/Positions.txt')
%hej=textscan(fid,['%*s %*s %*s ' repmat('%f %f %f',1,27) ' \n'],'headerlines',0);
hej=textscan(fid,['%*s %s %*s ' repmat('%f ',1,27*3)],'headerlines',0,'Whitespace',' \b\t\n\r','CollectOutput',true);
%core=textscan(sprintf,'VirtualSatelliteTraces%f.txt')
x=unique(hej{2}(:,1:3:end));
y=unique(hej{2}(:,2:3:end));
z=unique(hej{2}(:,3:3:end));
Ps=[reshape(hej{2}(:,1:3:end)',1,27*size(hej{2},1))',...
    reshape(hej{2}(:,2:3:end)',1,27*size(hej{2},1))',... 
    reshape(hej{2}(:,3:3:end)',1,27*size(hej{2},1))'];
for k=1:length(hej{1}) % collect cores
    core(k)=sscanf(char(hej{1}(k)),'VirtualSatelliteTraces%f.txt');
end
cores=reshape(repmat(core',1,27)',27*length(core),1);
%Pss=reshape(Ps',3,size(Ps,1)*size(Ps,2)/3)';
%%
fid2=fopen('/Users/Cecilia/TRED46/DATA/gridpoints.txt'); %add ,'a' to change file
fprintf(fid2,'%5g\n', [x;y;z]);
fclose(fid2);


%% Find core of positions between limits
xlims=[5.5 6.5];ylims=[8 9];zlims=[0 10];
xs=intersect(x(find(x>xlims(1))),x(find(x<xlims(2))));
ys=intersect(y(find(y>ylims(1))),y(find(y<ylims(2))));
zs=intersect(z(find(z>zlims(1))),z(find(z<zlims(2))));

% make combinations of vectors
combs=zeros(length(xs)*length(ys)*length(zs),3);
index=0;
for k=1:length(xs)
    for l=1:length(ys)
        for m=1:length(zs)
            index=index+1;
            combs(index,:)=[xs(k) ys(l) zs(m)];
        end
    end
end

[~,Psi,~]=intersect(Ps,combs,'rows');
sc=repmat(tocolumn(1:27),length(core),1);
%% Plot location of cores and satellite
figure;
axis([xlims(1) xlims(2) ylims(1) ylims(2) zlims(1) zlims(2)]);
for k=1:length(Psi)
    %plot3(Ps(Psi(k),1),Ps(Psi(k),2),Ps(Psi(k),3)); hold on;
    text(Ps(Psi(k),1),Ps(Psi(k),2),Ps(Psi(k),3),[num2str(cores(Psi(k))),'-',num2str(sc(Psi(k)))])
end

%% Plot E-fields of read data
for m=4
for k=1:3
    figure; h=irf_plot(9);
    for l=1:9
        irf_plot(h(l),[(1:17531)' squeeze(Ne(m,(k-1)*9+l,:,:))]);
        eval(['ylabel(h(l),''',num2str(core(m)),'/',num2str((k-1)*9+l),''')']);
    end
end
end

%% Calculate Epar
ncore=8;
Jepar=zeros(ncore,27,17531,1);
for l=1:ncore
    for k=1:27
        epar_time=irf_dot([(1:17531)' squeeze(Je(l,k,:,:))],...
                irf_norm([(1:17531)' squeeze(B(l,k,:,:))]));
            Jepar(l,k,:,:)=epar_time(:,2);
    end
end
        
%% Find distance along B-field
% match core 437:25-26
cr=4; sc=[11 12];
time=15830;
dist=torow(squeeze(Position(cr,sc(2),:)-Position(cr,sc(1),:)));
Bav=mean([torow(squeeze(B(cr,sc(2),time,:)));torow(squeeze(B(cr,sc(1),time,:)))]);
distB=dot(Bav,dist)/norm(Bav);
% Find timelag
[xcorrs, lags]=xcorr(Epar(2,26,:),Epar(2,25,:),50);
[maxcorr, lagind]=max(squeeze(xcorrs));
tlag=lags(lagind);
plot(lags,squeeze(xcorrs))
% Find velocity, distance/timelag
irf_units
velocity=distB/tlag;
vkms=velocity*Units.c*1e-3; %km/s

%% Create overview plot 
h=irf_plot(10);
isub=1;
t=tocolumn(squeeze(1:size(E,3)));
sc(1)=12;

if 1 % B
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(B(cr,sc(1),:,:))]);
    ylabel(hca,'B')
end
if 1 % E
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(E(cr,sc(1),:,:))]);
    ylabel(hca,'E')    
end
if 1 % Epar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Epar(cr,sc(1),:))]);
    ylabel(hca,'E_{||}')    
end
if 1 % Ji
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ji(cr,sc(1),:,:))]);
    ylabel(hca,'J_i')
end
if 1 % Je
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Je(cr,sc(1),:,:))]);    
    ylabel(hca,'J_e')
end
if 1 % Ni
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ni(cr,sc(1),:))]);
    ylabel(hca,'N_i')
end
if 1 % Ne
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ne(cr,sc(1),:))]);
    ylabel(hca,'N_e')
end
if 1 % Ntot
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ne(cr,sc(1),:))+squeeze(Ni(cr,sc(1),:))]);
    ylabel(hca,'N_{tot}')
end
if 1 % Ve
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t Units.c*1e-3*squeeze(Je(cr,sc(1),:,:))./repmat(squeeze(Ne(cr,sc(1),:)),1,3)]);
    ylabel(hca,'V_e [km/s]')
end
if 1 % Vi
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t Units.c*1e-3*squeeze(Ji(cr,sc(1),:,:))./repmat(squeeze(Ni(cr,sc(1),:)),1,3)]);
    ylabel(hca,'V_i [km/s]')
end      
%% Read data and plot location of eh observations and |B|/B
fid=fopen('/Users/Cecilia/TRED46/DATA/ehpoints.txt')
data=textscan(fid,'%f %f %*[^\n]','delimiter','\n');
figure;axis([5 8 7 10 0 10]);
xlabel('x (x_{GSM})');ylabel('y (z_{GSM} sign?)');zlabel('z (-y_{GSM} sign?)');
for k=1:length(data{1})        
    cor=find(core==data{1}(k));
    sat=data{2}(k);    
    plot3(Bnorm(1),Bnorm(2),Bnorm(3)); hold on;
    text(Position(cor,sat,1),Position(cor,sat,2),Position(cor,sat,3),[num2str(data{1}(k)),'-',num2str(data{2}(k))])
end
%% Integrate E
Phipar=zeros(size(Epar));
for k=1:length(core)
    for p=1:27
        phi=irf_integrate([t squeeze(Epar(cr,sc(1),:,:))]);
        Phipar(k,p,:,:)=phi(:,2);
    end
end
%% Create match plot 
figure;
h=irf_plot(4);
isub=1;
t=tocolumn(squeeze(1:size(E,3)));
if 1
    hca=h(isub)
    isub=isub+1;
    irf_plot(hca,[t squeeze(B(cr,sc(1),:,:))]);
    ylabel(hca,'B')
end
if 1
    hca=h(isub) 
    isub=isub+1;
    irf_plot(hca,{[t squeeze(Epar(cr,sc(1),:,:))], [t squeeze(Epar(cr,sc(2),:,:))], [t squeeze(Epar(cr,sc(2),:,:))]},'comp','dt',[0 0 1*tlag]);
    ylabel(hca,'E_{||}')    
end
if 1
    hca=h(isub) 
    isub=isub+1;
    irf_plot(hca,{[t squeeze(Phipar(cr,sc(1),:,:))], [t squeeze(Phipar(cr,sc(2),:,:))], [t squeeze(Phipar(cr,sc(2),:,:))]},'comp','dt',[0 0 1*tlag]);
    ylabel(hca,'Phi_{||}')    
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ji(cr,sc(1),:,:))]);
    ylabel(hca,'J_i')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Je(cr,sc(1),:,:))]);    
    ylabel(hca,'J_e')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ni(cr,sc(1),:))]);
    ylabel(hca,'N_i')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ne(cr,sc(1),:))]);
    ylabel(hca,'N_e')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t squeeze(Ne(cr,sc(1),:))+squeeze(Ni(cr,sc(1),:))]);
    ylabel(hca,'N_{tot}')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t Units.c*1e-3*squeeze(Je(cr,sc(1),:,:))./repmat(squeeze(Ne(cr,sc(1),:)),1,3)]);
    ylabel(hca,'V_e [km/s]')
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t Units.c*1e-3*squeeze(Ji(cr,sc(1),:,:))./repmat(squeeze(Ni(cr,sc(1),:)),1,3)]);
    ylabel(hca,'V_i [km/s]')
end  

%% Find velocity with c_4_v
figure;
h=irf_plot(4);
t=tocolumn(squeeze(1:size(E,3)));
cr=444*[1 1 1 1];
sc=[10 11 12 13];
for k=1:4;
    hca=h(k);
    irf_plot(hca,[t squeeze(E(cr(k),sc(k),:,:))]);
    ylabel(hca,'E')
end




