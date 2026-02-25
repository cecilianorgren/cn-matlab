function ln = ln2(t1,t2,gsmB3,gsmB4,dx,dy,M)
mu0=4*pi*1e-7;
c_eval('gseB?=c_coord_trans(''gsm'',''gse'',gsmB?,''CL_ID'',?);',3:4);

dB=irf_add(1,gseB3,-1,gseB4);
dB=cn_toepoch(t1,t2,dB);
gseB3=cn_toepoch(t1,t2,gseB3);
gseB4=cn_toepoch(t1,t2,gseB4);

prel=M*dB(:,2:4)';
bdB=[dB(:,1) prel'];
prel=M*gseB3(:,2:4)';
bB3=[gseB3(:,1) prel'];
prel=M*gseB4(:,2:4)';
bB4=[gseB4(:,1) prel'];

bB3=irf_resamp(bB3,bdB);
bB4=irf_resamp(bB4,bdB);

jk=[bdB(:,1) bdB(:,2)*dx*1e3*1e-9/mu0]; % A/m^2
jnorm=[bdB(:,1) bdB(:,3)*dy*1e3*1e-9/mu0];

jxB=irf_multiply(1,jk,1,bB3(:,[1 4]),1);


ln=[jxB(:,1) (0.1*30)./jxB(:,2)];
ln=[jxB(:,1) jxB(:,2)/(0.1*30)];
figure;irf_plot(jk);hold on;
irf_plot(bB3(:,[1 4]),'g');hold on;