% cold ions 
% cold beam electrons
% hot background electrons
disprel = @(opi,ope,opeb,vti,vte,vteb,om,k) (1 - opi^2./om^2 - opeb^2./((om-k*vteb).^2)-1i*sqrt(pi)*(om*ope./(k.^3*vte.^3))*exp(-(om./k/vte).^2));

%%

imo=@(reo,opekv)(sqrt(3*reo.^2-2*reo.*opekv));
zero=@(reo,opekv)(reo.^3-opekv.*reo.^2+opekv.^2.*reo/4-1/16);
%% Three beams, one ion, two opposite electron
myfun=@(o,k,opi,op1,op2,vi,v1,v2)(1+...
    -opi^2./(o-k.*vi).^2+...
    -op1^2./(o-k.*v1).^2+...
    -op2^2./(o-k.*v2).^2);
%% Find zeros
k=-2:0.1:1.5;
for l=1:length(k)
    imre(l)=fzero(@(o)myfun(o,k(l),1,1,1,1,1,1),0);
end

%%
myfun2=@(o,k,opi,op1,op2,vi,v1,v2)(1+opi^2./o^2+...
    +(o-k*v1).*(o-k*v2)./sqrt(op1^2*(o-k*v1)+op2^2*(o-k*v2)));
fzero(@(o)imag(myfun2(o,1,1,1,1,1,1,10)),0.1*i)

%%
myfun3=@(o,opekv)(o.^3-opekv.*o.^2+1/2);
fzero(@(o)imag(myfun3(o,0)),0)
%%
opekv=-2:0.05:1.5;
for l=1:length(opekv)
    reom(l)=fzero(@(o)zero(o,opekv(l)),0);
end
imom=imo(reom,opekv);
plot(opekv,reom,opekv,imom)
xlabel('kv_d-\omega_{pe}','fontsize',14)
ylabel('\omega_r , \gamma   [\omega_{pi}^{2/3}\omega_{pe}^{1/3}]','fontsize',14)
title('Buneman instability, cold plasma limit','fontsize',14)

%% Compare for slight spread in kvd-ope
% kvd=ope
one=@(or)(8*or.^3-1/2);
twoplus=@(or)(8*or.^3-or.^2-1/2);
twominus=@(or)(8*or.^3+or.^2-1/2);

or=-5:0.1:5;

plot(or,one(or),or,twoplus(or),or,twominus(or))

fzero(one,0)
fzero(twominus,0)
fzero(twoplus,0)
%% Phase velocity
imo2=@(reo,opekv)(sqrt(3*reo.^2-2*reo.*opekv));
zero2=@(reo,opekv)(reo.^3-opekv.*reo.^2+opekv.^2.*reo/4-1/16);
vph=@(spread)(fzero(@(o)zero(o,spread),0)*(Units.me/Units.mp)^(2/3)/(1-(Units.me/Units.mp)^(2/3)*spread^(1/3)));
%% Phase velocity dependant on spread.
vph=@(spread)(fzero(@(o)zero(o,spread),0).*(Units.me/Units.mp)^(2/3)./(1+spread.*(Units.me/Units.mp)^(2/3)));
spread=-2:0.1:1.5;
for k=1:length(spread)
    vph_vec(k)=vph(spread(k));
end
plot(spread,vph_vec)

%% Find for fractional R electron drift
funct=@(R,ope,opi,kvd,om)(R*ope^2*om.^2+(opi^2+ope^2*(1-R)).*(kvd-om).^2-om.^2.*(kvd-om).^2);
fpe=3000; % Hz
ope=2*pi*fpe; % Rad/s
fpi=300; % Hz
opi=2*pi*fpi; % Rad/s
kvd=1e2:1e1:2e3;
R=0.0;%0.5;

for k=1:length(kvd)
    omega(k)=fzero(@(om)funct(R,ope,opi,kvd(k),om),1000);
end

plot(kvd,omega)
clear omega

%% Two beams with background
imo=@(reo,opekv)(sqrt(3*reo.^2-2*reo.*opekv));
zero=@(reo,opekv)(reo.^3-opekv.*reo.^2+opekv.^2.*reo/4-1/16);