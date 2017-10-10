function cn_cjeck_j(tint,gsmB3,gsmB4,gsePos3,gsePos4,flag)
% %%%%% check j %%%%%%
% inputs: tint, gsmB3, gsmB4, gsePos3, gsePos4, flag
% flag = 1 - fixed coordinate system during the whole time interval
% flag = -1 - the coordinate system moves with b for each time step

t1=tint(1); t2=tint(2);

% b-field
b3=irf_tlim(gsmB3(:,1:4),t1,t2);   % b3
b4=irf_tlim(gsmB4(:,1:4),t1,t2);   % b4
b=irf_add(0.5,b3,0.5,b4);   % bav
gsmdb=irf_add(1,b3,-1,b4);  % db

if flag % fixed coordinate system during the whole time interval
    
else % the coordinate system moves with b for each time step
    
end

%% Clear variables
%clear b3 b4 b gsmdb z zav sagsm3 sagsm4 sagsmav spgsmav facb3 facb4 facdb
%clear gsmdist facdist facdistav facj gsmj
%%
%t1=toepoch([2007 08 31 10 13 10]); t2=toepoch([2007 08 31 10 13 50]);
%t1=toepoch([2007 08 31 10 19 00]); t2=toepoch([2007 08 31 10 19 10]);
%%


z=irf_norm(b(:,1:4));       % zav for each time step
zav=mean(z(2:4),1);         % zav during time interval                        
c_eval('sagsm?=c_coord_trans(''dsi'',''gsm'',[z(:,1) repmat([0 0 1],size(b?,1),1)/sqrt(2)],''cl_id'',?);',3:4)
                            % spin axis for each sc and each time step
sagsm=[sagsm3(:,1) (sagsm3(:,2:4)+sagsm4(:,2:4))/2]; % spin axis average between sc
sagsmav=mean(sagsm(:,2:4)); % spin axis average during time interval 
spgsmav=irf_cross(sagsmav,zav); % average vector perp to b and in spin plane (spav)
spgsm=irf_cross(sagsm,z);
facb3=irf_lmn(b3,zav,spgsmav,'L');
facb4=irf_lmn(b4,zav,spgsmav,'L');
facdbav=irf_lmn(gsmdb,zav,spgsmav,'L'); % fac delta b based on bav and spav
facdb=irf_lmn(gsmdb,z,spgsmav,'L'); % fac delta b different for each time step
facb3=irf_abs(facb3);
facb4=irf_abs(facb4);
%% positions
c_eval('gsmPos?=c_coord_trans(''gse'',''gsm'',cn_toepoch(t1,t2,gsePos?),''cl_id'',?);',3:4)
gsmdist=irf_add(1,gsmPos3,-1,gsmPos4);
facdistav=irf_lmn(gsmdist,zav,spgsmav,'L'); % gives back x (B) y (third) z (spin plane)
facdist=irf_lmn(gsmdist,z,spgsm,'L');

%% current
facj=cn_j(facdb,facdist); % A
gsmj=cn_j(gsmdb,gsmdist); % A
facj(:,2:4)=facj(:,2:4)*1e9;
gsmj(:,2:4)=gsmj(:,2:4)*1e9;

%% plot
figure(33);h=irf_plot(7);isub=1;
hca=h(isub);isub=isub+1;
irf_plot(hca,{facb3(:,[1 2]),facb4(:,[1 2])},'comp');
ylabel(hca,'B_x  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,{facb3(:,[1 3]),facb4(:,[1 3])},'comp');
ylabel(hca,'B_y  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,{facb3(:,[1 4]),facb4(:,[1 4])},'comp');
ylabel(hca,'B_z  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,{facb3(:,[1 5]),facb4(:,[1 5])},'comp');
ylabel(hca,'|B|  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,facdb(:,1:4));ylabel(hca,'\Delta B  [nT]');irf_legend(hca,{'x (B)','y','z (spin plane)'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,facdistav);ylabel(hca,'s/c separation  [km]');irf_legend(hca,{'x (B)','y','z (spin plane)'},[0.02 0.04])
hca=h(isub);isub=isub+1;
irf_plot(hca,facj);ylabel(hca,'j  [nA]');irf_legend(hca,{'x (B)','y','z (spin plane)'},[0.02 0.04])
irf_zoom(h,'x',[facj(1,1) facj(end,1)]);
%% print
if 1 % make date/time strings
    t1str_p=datestr(epoch2date(t1),'HHMMSSFFF');
    t2str_p=datestr(epoch2date(t2),'MMSSFFF');
end
eval(['print -depsc2 /Users/Cecilia/Dropbox/Cecilia/EH/current/',t1str_p,'-',t2str_p,'_j.eps']);
