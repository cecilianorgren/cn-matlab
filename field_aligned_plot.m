%% Load interval
tint=[toepoch([2007 08 31 09 45 00]) toepoch([2007 08 31 12 15 00])];

%% Take out B for interval with two biggest blops
c_eval('[caagseB?,~,gseB?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',3:4);
tint_zoom=[toepoch([2007 08 31 10 18 42]) toepoch([2007 08 31 10 18 44])];
i_start=find(gseB3(:,1)>tint_zoom(1),1,'first');
i_end=find(gseB3(:,1)<tint_zoom(2),1,'last');
c_eval('gseB?zoom=gseB?(i_start:i_end,:);',3:4);

%% Even smaller interval
% Nint=size(gseB3zoom,1);
% c_eval('gseB?zooom=gseB?zoom(fix(4*Nint/10):fix(Nint/2),:);',3:4);

%% Take average value of B
N=size(gseB3zoom,1);
c_eval('gseB?av=[sum(gseB?zoom(:,2))/N sum(gseB?zoom(:,3))/N sum(gseB?zoom(:,4))/N];',3:4);

%% Defining B hat
c_eval('B?mag = mag(gseB?av);',3:4);
c_eval('gseB?hat = gseB?av/B?mag;',3:4);

%% Taking out xhat and yhat
c_eval('z?hat=gseB?hat;',3:4);
c_eval('x?hat=cross(cross(z?hat,[0 1 0]),z?hat);',3:4);
c_eval('x?hat=x?hat/mag(x?hat);',3:4);
c_eval('y?hat=cross(z?hat,x?hat);',3:4);

%% plot
plot3([0 x3hat(1)],[0 x3hat(2)],[0 x3hat(3)],'--b'); hold on;
plot3([0 y3hat(1)],[0 y3hat(2)],[0 y3hat(3)],'--g'); hold on;
plot3([0 z3hat(1)],[0 z3hat(2)],[0 z3hat(3)],'--r'); hold on;
legend('y_{B}','x_{B}','B')

plot3([0 1],[0 0],[0 0],'b'); hold on;
plot3([0 0],[0 1],[0 0],'g'); hold on;
plot3([0 0],[0 0],[0 1],'r'); hold on;
xlabel('x_{gse}');ylabel('y_{gse}');zlabel('z_{gse}');

%% Transformation matrix
M=[x3hat;y3hat;z3hat];

%% Transforming from gse to field aligned coordinate system 
% xnew=M*xold'

%% Plotting position of C3 C4 in field aligned 2D system
c_eval('[caaPos?,~,gsePos?]=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'');',3:4);

%bcPos=[gsePos(:,1) M*gsePos(1,2:4)];
c_eval('gsePos?=gsePos?(fix((i_end+i_start)/2),:);',3:4);
%%
c_eval('gsePos?=gsePos?;',3:4);


%% %% Plotting spacecraft trajectory

x1=gsePos(1,2:4)';
x2=gsePos(fix(end/4),2:4)';
x3=gsePos(fix(end/2),2:4)';
x4=gsePos(fix(3*end/4),2:4)';
x5=gsePos(end,2:4)';

x0=(x1+x2+x3+x4+x5)/5;
for n=1:5
    eval(['x',num2str(n),'=x',num2str(n),'-x0;'])
    eval(['x',num2str(n),'b=M*x',num2str(n),';'])
end

x0b=M*x0;

%plot3([x1(1) x2(1) x3(1) x4(1) x5(1)],[x1(2) x2(2) x3(2) x4(2) x5(2)],[x1(3) x2(3) x3(3) x4(2) x5(2)])
plot3([x1b(1) x2b(1) x3b(1) x4b(1) x5b(1)],[x1b(2) x2b(2) x3b(2) x4b(2) x5b(2)],[x1b(3) x2b(3) x3b(3) x4b(2) x5b(2)],'b')
hold on;
xlabel('x_{B} [km]');ylabel('y_{B} [km]');zlabel('B_{||} [km]');

%% Plotting ion velocities
c_eval('[caaVi?,~,Vi?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);
%c_eval('dionv?=c_coord_trans(''gse'',''dsi'',ionv?,''CL_ID'',4);',3:4);   
        





