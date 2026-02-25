%%
clear lambda kvec chii chie1 chie2
lambda=(5:5000)'*1e3; % m
kvec=2*pi./lambda; % m^-1

R=1;%0.0015; % electron beam fraction in density

vi=1000e3;
ve1=20000e3;
ve2=1000e3;

vti=700e3;
vte1=15000e3;
vte2=20000e3;

opi=400; % rad/s
ope=5e3; % rad/s
ope1=R*ope; % rad/s
ope2=(1-R)*ope; % rad/s

omega=linspace(ope,opi,10);

for p=1:length(kvec)    
    for q=1:length(omega)
        k=kvec(p);
        etai=(omega(q)-k.*vi)./k./vti;
        etae1=(omega(q)-k.*ve1)./k./vte1;
        etae2=(omega(q)-k.*ve2)./k./vte2;

        chii(p,q)=(1+etai*cef(etai,16))*2*opi.^2./k./vti.^2;
        chie1(p,q)=(1+etae1*cef(etae1,16))*2.*ope1.^2./k./vte1.^2;
        chie2(p,q)=(1+etae2*cef(etae2,16))*2.*ope2.^2./k./vte2.^2;
    end
end

one=ones(length(kvec),length(omega));
eps=one+chii+chie1+chie2;

plot(kvec,real(eps),kvec,imag(eps))
%figure(1);surf(kvec,omega,real(eps))
%figure(2);surf(kvec,omega,imag(eps))

%% Test to find zeros
f=@(x)(x.^2-2*i);
fsolve(f,2+2*i)

%% Cold plasma buneman
buneman=@(o,k,opi,ope,R,vi,ve1,ve2)...
    (1+opi^2./(o-k*vi).^2+R*ope^2./(o-k*ve1)+(1-R)*ope^2./(o-k*ve2));

nk=50;    
R=1;
vi=0;
ope=2.20e3*2*pi; % rad/s
%k=ope/ve1;
opi=ope*sqrt(1/1836); % rad/s
ve2=0;
ve1=10000e3; %m/s

k=linspace(ope-ope^(1/3)*opi^(2/3)*2,ope+ope^(1/3)*opi^(2/3)*1.5,nk)/ve1;

os=ope^(1/3)*opi^(2/3);
ks=ope/ve1;

%k=linspace(opi,ope,nk);
for p=1:nk
    omega(p)=fsolve(@(o)(buneman(o,k(p),opi,ope,R,vi,ve1,ve2)),ope^(1/3)*opi^(2/3));
    gamma(p)=fsolve(@(o)(buneman(o,k(p),opi,ope,R,vi,ve1,ve2)),1*ope^(1/3)*opi^(2/3));
end
plot(k/ks,omega/os,k/ks,gamma/os)