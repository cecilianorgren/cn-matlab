% despin 42 (24?)
c_load('wE3p34')
c_load('wE3p32')
c_load('Atwo3')
dwE3p24 = irf_add(1,wE3p34,-1,wE3p32);
es = [dwE3p24 dwE3p24(:,1)*0];
coef=[[1 0 0];[1 0 0]];
phase = c_phase(es(:,1),Atwo3);

e = c_efw_despin(es,phase,coef,'asym42');