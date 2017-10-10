%heeeeeeej
[x y z correlation PhiE Bz Ek En ufEn ufEk] = tool_direction(cn_toepoch(t1,t2,gsmB3),cn_toepoch(t1,t2,gsmE3));
index=find(correlation==max(correlation));
%% 
n=0.1;
[v correlation_amp] = tool_velocity(15,Bz(:,2),PhiE(:,1+index),n,100);
index_v=find(correlation_amp==max(correlation_amp));
v(index_v)
