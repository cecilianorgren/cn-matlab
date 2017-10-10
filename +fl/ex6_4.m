% ex 6.4

%% vz(r)
R1=1;
R2=6;
dpdz=-1.5;
eta=1;
vz = @(r)(-(R1^2*dpdz/4/eta)*(1-r.^2/R1/R1-(1-R2^2/R1/R1).*log(r/R1)/log(R2/R1)));
r=linspace(R1,R2,50);
plot(vz(r),r,vz(r),r)
title(['R1 = ' num2str(R1) ', R2 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
ylabel('v_z(r)')
xlabel('r')

%% Q
R1=1;
R2=linspace(R1,3*R1,50);
Q = @(R1,R2)(pi*dpdz*(R1^2-R2.^2)/8/eta.*(R1^2+R2.^2+(R1^2-R2.^2)./log(R2./R1)));
area=pi*(R2.^2-R1.^2);
plot(R2,Q(R1,R2))
title(['R1 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
ylabel('Q(R2)')
xlabel('R2')

%%
plot(R2,Q(R1,R2)./area)
title(['R1 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
ylabel('Q(R2)/area')
xlabel('R2')
%% v_max
nr=100;
R1=1;
nr2=50;
R2=linspace(R1,5*R1,nr2);
vmax=zeros(nr2,1);
vz = @(r,R1,R2)(-(R1^2*dpdz/4/eta)*(1-r.^2/R1/R1-(1-R2.^2/R1/R1).*log(r/R1)./log(R2/R1)));
for n=1:nr2;    
    r=linspace(R1,R2(n),nr);
    vmax(n)=max(vz(r,R1,R2(n)));
end

plot(R2,vmax)
title(['R_1 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
ylabel('v_{max}')
xlabel('R_2')
%vmax=vz((R1+R2)/2,R1,R2);
%% v_mean
nr=100;
R1=1;
nr2=50;
R2=linspace(R1,5*R1,nr2);
vmean=zeros(nr2,1);
vz = @(r,R1,R2)(-(R1^2*dpdz/4/eta)*(1-r.^2/R1/R1-(1-R2^2/R1/R1).*log(r/R1)/log(R2/R1)));
for n=1:nr2;    
    r=linspace(R1,R2(n),nr);
    vmean(n)=mean(vz(r,R1,R2(n)));
end
plot(R2,vmean,R2,vmax)
title(['R_1 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
ylabel('v')
xlabel('R_2')
legend({'v_{mean}','v_{max}'},'location','northwest')

%%

%% vz(r)
R1=1;
R2=2;
dpdz=-1.5;
eta=1;

% v(r)
hca = subplot(2,1,1);
vz = @(r,R1,R2)(-(R1^2*dpdz/4/eta)*(1-r.^2/R1/R1-(1-R2^2/R1/R1).*log(r/R1)/log(R2/R1)));
R2a = 2;
R2b = 4;
ra = linspace(R1,R2a,50);
rb = linspace(R1,R2b,50);

plot(hca,vz(ra,R1,R2a),ra/R1,vz(rb,R1,R2b),rb/R1)
%title(['R_1 = ' num2str(R1) ', R_2 = ' num2str(R2) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
title('Velocity profile')
xlabel('v_z(r)')
ylabel('r/R_1')
legend('R_2=2R_1','R_2=4R_1')%,'location','best')
hca.FontSize = 14;

% Q(R2)
hca = subplot(2,1,2);
R1=1;
R2=R1*linspace(1,3,50);
Q = @(R1,R2)(pi*dpdz*(R1^2-R2.^2)/8/eta.*(R1^2+R2.^2+(R1^2-R2.^2)./log(R2./R1)));
area=pi*(R2.^2-R1.^2);
plot(hca,R2/R1,Q(R1,R2))
%title(['R_1 = ' num2str(R1) ', dp/pz = ' num2str(dpdz) ', \eta = ' num2str(eta)])
title('Volume flow rate')
ylabel('Q(R_2)')
xlabel('R_2/R_1')
hca.FontSize = 14;