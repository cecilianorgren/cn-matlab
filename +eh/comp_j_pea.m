% comparing current to electron distributions
cd /Users/Cecilia/Data/BM/20070831
tint = toepoch([2007 08 31 10 17 10;...
                2007 08 31 10 18 10])';
%% Magnetic field
load matlabdiEB ;
%% Spacecraft positions
c_eval('gsePos?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);         
c_eval('diPos?=c_coord_trans(''gse'',''isr2'',gsePos?,''cl_id'',?);',3:4); 
%% Make field aligned coordinate system
z = irf_norm(cn.mean([irf_tlim(diB3(:,1:4),tint); irf_tlim(diB4(:,1:4),tint)],1));
disp(['z=' cn.vector_to_string(z)])
y = irf_norm(cross(z,[2 0 0]));
x = cross(y,z);
c_eval('facB? = irf_newxyz(diB?,x,y,z);',3:4);
c_eval('facPos? = irf_newxyz(diPos?,x,y,z);',3:4);
% Calculate current
facj = cn_j(irf_add(-1,facB3,1,facB4),irf_add(-1,facPos3,1,facPos4)); % A
% Plot current
irf_plot(irf_tappl(facj,'*1e9')) % nA
irf_zoom('x',tint)
irf_legend({'x_{fac}','y_{fac}','z_{fac}'},[0.05 0.9])
ylabel('j [nA]')
title(['x=' cn.vector_to_string(x) '  y=' cn.vector_to_string(y) '  z=' cn.vector_to_string(z)])

 %   c_eval('cn.print(''current_?'',''current'');',ii)
end

%% Many
for ii=1:8
    ax{ii}=subplot(4,2,ii);

z = irf_norm(cn.mean([irf_tlim(diB3(:,1:4),tint); irf_tlim(diB4(:,1:4),tint)],1));
y = irf_norm(cross(z,[4 ii 0]));
x = cross(y,z);
c_eval('facB? = irf_newxyz(diB?,x,y,z);',3:4);
c_eval('facPos? = irf_newxyz(diPos?,x,y,z);',3:4);
% Calculate current
facj = cn_j(irf_add(-1,facB3,1,facB4),irf_add(-1,facPos3,1,facPos4)); % A
% Plot current
irf_plot(ax{ii},irf_tappl(facj,'*1e9')) % nA
irf_zoom(ax{ii},'x',tint)
irf_legend(ax{ii},{'x_{fac}','y_{fac}','z_{fac}'},[0.05 0.9])
ylabel(ax{ii},'j [nA]')
title(ax{ii},['x=' cn.vector_to_string(x) '  y=' cn.vector_to_string(y) '  z=' cn.vector_to_string(z)])
% Mark when HEEA is into B

 %   c_eval('cn.print(''current_?'',''current'');',ii)
end