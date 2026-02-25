c_eval('gsex?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 1 0 0],''CL_ID'',?);',3:4);
c_eval('gsey?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 0 1) 0],''CL_ID'',?);',3:4);

c_eval('gseEx?=c_coord_trans(''dsi'',''gse'',[Ee3(1:2) 0 0],''CL_ID'',?);',3:4);
c_eval('gseEy?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 0 Ee3(3) 0],''CL_ID'',?);',3:4);

%%
figure
s=30;
m=-110;
plot3(s*[0 1],s*[0 0],s*[0 0],'--c'); hold on;
plot3(s*[0 0],s*[0 1],s*[0 0],'--c'); hold on;
plot3(s*[0 0],s*[0 0],s*[0 1],'--c'); hold on;


plot3([0 0],[0 11.6],[0 -0.3],'b'); hold on;
plot3([0 52.8],[0 -0.2],[0 -6.8],'g'); hold on;
plot3(m*[0 -0.8],m*[0 -0.1],m*[0 -0.6],'r'); hold on;

%%
cn_scalar(Bh,gsex3(2:4))
cn_mag(gseEx3(2:4))
