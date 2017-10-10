% local expansion
ve = cn_eV2v(1000,'eV');
v = linspace(cn_eV2v(0.1,'eV'),cn_eV2v(5e3,'eV'),1000);
f = cn.maxwellian(v,35,20,0,'e','3D'); % f = cn.maxwellian(v,T,n,vd,species,optional);      
fd = cn.maxwellian(v,35,20,cn_eV2v(30,'eV'),'e','3D'); 
en = -100;
oce = 1;
flexp = 1*(1-en*v'/max(v)).*cn.maxwellian(v,35,20,0,'e','3D'); 
semilogx(v,f,v,fd,v,flexp)
   
