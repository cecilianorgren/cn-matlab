function [dt dx dy] = cn_tool_dt(E3,E4,v,t1,t2,M,gsePos3_raw,gsePos4_raw)

c_eval('E?=cn_toepoch(t1,t2,E?);',3:4)

c_eval('gsePos?=cn_toepoch(t1,gsePos?_raw);',3:4)
gsePos3(2:4);
M(2,:);
c_eval('gsePos?y=M(2,:)*gsePos?(2:4)'';',3:4)
c_eval('gsePos?x=M(1,:)*gsePos?(2:4)'';',3:4)
dx=gsePos4x-gsePos3x;
dy=gsePos4y-gsePos3y;

dt=dy/v;

% Transforming to field aligned system        
Pos3=M*gsePos3(2:4)';
Pos4=M*gsePos4(2:4)';
Pos0=(Pos3+Pos4)/2;
c_eval('Pos?c=Pos?-Pos0;',3:4);

boxside=max(abs([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]));
plot(Pos3c(1),Pos3c(2),'go','markersize',13,'linewidth',1.5); hold on;
plot(Pos4c(1),Pos4c(2),'bv','markersize',13,'linewidth',1.5); hold on;
c3=text(Pos3c(1),Pos3c(2),'    C3','color','k');
c4=text(Pos4c(1),Pos4c(2),'    C4','color','k');

normdist=abs(Pos4c(1)-Pos3c(1));
zdist=Pos4c(3)-Pos3c(3);
kdist=abs(Pos4c(2)-Pos3c(2))
cn_mag(Pos3c-Pos4c);

xlabel('x   [km]'); ylabel('y   [km]');       
axis([-boxside boxside -boxside boxside]*1.1)
axis square       

%cn_plot_corr(E3,E4,dt,'tool');