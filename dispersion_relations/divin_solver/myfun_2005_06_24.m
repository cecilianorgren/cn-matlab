function [ output_args ] = myfun_dispersion5A( input_args,kk,fv)
% Attention: takes a complex argument as a (real(x), imag(x)) vector
% A special version, use Yoon/Lui (2008) parameters to study Buneman-LHDI
% interaction.

% fv = [omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2,lamD];
omega_pi = fv(1);
omega_pe1 = fv(2);
omega_pe2 = fv(3);
vthi = fv(4);
vthe1 = fv(5);
vthe2 = fv(6);
vdi = fv(7);
vde2 = fv(8);
vde2 = fv(9);
lamD = fv(10);

%%
z=input_args(1)+1i*input_args(2);


kx = kk(1)*lamD;  %omega_ci/vthi;%omega_ce/vthe;
ky = kk(2)*lamD;  %*omega_ci/vthi;%omega_ce/vthe;
kz = kk(3)*lamD;%  - (1.2*solve_parameterZ*omega_ce/vthe); %*omega_ci/vthi;%*omega_ce/vthe;
kk=sqrt(kx.^2 + ky.^2 - kz.^2);


%qq=solve_parameter*sqrt(mm);

bb=((kx*sin(alpha_angle))-(ky*cos(alpha_angle)) )*vde*cos(alpha_angle)/omega_ce;
lam=((  vthe*vthe*  ( -kz.^2 + ( (kx*sin(alpha_angle)) - ( ky*cos(alpha_angle))   ).^2) )/(2*omega_ce*omega_ce)) + ((kz*vde*cos(alpha_angle))/omega_ce) ;
chi=(z/vthe/( kx*cos(alpha_angle)  +  ky*sin(alpha_angle)    )  )  -  ( vde*sin(alpha_angle) / vthe );

return_residual= 1 - ...
    ((omega_pi*omega_pi/(kk*kk*vthi*vthi))*Zplasma_d(( z - ky*vdi)/(kk*vthi)  ))+...
    (2*omega_pe^2/(kk*kk*vthe*vthe) ) * (...
1+...
((z - ky*vde)/(( kx*cos(alpha_angle) +  ky*sin(alpha_angle) )*vthe)*...
besselj(0,bb)*besseli(0,lam)*exp(-lam)*Zplasma(chi) ...
)...
);

if (isinf(return_residual)) return_residual=1000;end;
if (isnan(return_residual)) return_residual=1000;end;

output_args=[real(return_residual), imag(return_residual)];

end

