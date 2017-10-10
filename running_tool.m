c_eval('gsmE?=c_coord_trans(''gse'',''gsm'',gseE?,''cl_id'',?);',3:4)
%%
c_eval('gsmPos?=c_coord_trans(''gse'',''gsm'',gsePos?,''cl_id'',?);',3:4)
%%
vvv=-1500;

[xc yc zc corrx corramp]=cn_tool_v(gsmB3,gsmE3,N*1e6,N*1e6,gsmPos4,vvv,t1,t2,100);

ix=find(corrx==max(corrx));
Mx=[xc(ix,:); yc(ix,:); zc(ix,:)];
propx=xc(ix,:);
normx=yc(ix,:);

iamp=find(corramp==min(corramp));
Mamp=[xc(iamp,:); yc(iamp,:); zc(iamp,:)];
propamp=xc(iamp,:);
normamp=yc(ix,:);
%%
figure;
[ttt xxx yyy]=cn_tool_dt(gsmE3,gsmE4,1300,t1,t2,Mamp,gsmPos3,gsmPos4)
%xchoice=x(6,:);

%%
%d=cn_event_conf(t1,t2,gseB3,gseB4,codifVi3,codifVi4,gseExB3,gseExB4,gsePos3,gsePos4,eVTe3,eVTi4,peaNeav,Niav,scpNe3,peaNe3,xchoice,t1str,t2str,'poster');

t0choice=ttt;
vchoice=vvv;
%%
cn_plot_corr(cn_toepoch(t1,t2,gseE3),cn_toepoch(t1,t2,gseE4),-ttt,'mycorr');
