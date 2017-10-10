function [ output_args ] = myfun_dispersion5A( input_args,solve_parameterX,solve_parameterY,solve_parameterZ )
% Attention: takes a complex argument as a (real(x), imag(x)) vector
% A special version, use Yoon/Lui (2008) parameters to study Buneman-LHDI
% interaction.

% fv = [omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2,lamD];
omega_pi = 

%% execute that in Matlab before running solver script.
mm=1836;                       % ion to electron mass ratio (1836 default)
alpha_angle=0.00001;           % need to be a small value (no guide field, see Yoon/Lui 2008 paper
 epsilon = 1.21*10^(-5);       % epsilon = (omega_ci/omega_pi)^2

%epsilon=10^(-8);

%omega_pi = 658.3517;
 omega_pi = 721;               % ion angular plasma frequency, in seconds^{-1}
 
 vthi = 600000;                % ion thermal velocity, m/s
%vthi= 875400;
 
% 180 000
vdi = 180000;                  % ion drift velocity, m/s
omega_pe = omega_pi*sqrt(mm);  % electron angular plasma frequency
 vthe=20000000;                % electron thermal velocity, m/s
%%%% vthe = vthi*sqrt(ttt*mm);
 
 vde = -150000;                % electron drift velocity, m/s
 omega_ci=omega_pi*sqrt(epsilon);
 omega_ce=mm*omega_ci;
 omega_lh=sqrt(mm)*omega_ci;
%%
 
 


 
 

z=input_args(1)+1i*input_args(2);


kx = solve_parameterX*omega_ce/vthe;  %omega_ci/vthi;%omega_ce/vthe;
%ky=solve_parameter*omega_lh/vthi;
ky = solve_parameterY*omega_ce/vthe;  %*omega_ci/vthi;%omega_ce/vthe;
%kz=0.5*ky*1i;
kz = -1i*solve_parameterZ*omega_ce/vthe;%  - (1.2*solve_parameterZ*omega_ce/vthe); %*omega_ci/vthi;%*omega_ce/vthe;
kk=sqrt(kx.^2 + ky.^2 - kz.^2);


%qq=solve_parameter*sqrt(mm);

bb=((kx*sin(alpha_angle))-(ky*cos(alpha_angle)) )*vde*cos(alpha_angle)/omega_ce;

lam=((  vthe*vthe*  ( -kz.^2 + ( (kx*sin(alpha_angle)) - ( ky*cos(alpha_angle))   ).^2) )/(2*omega_ce*omega_ce)) + ((kz*vde*cos(alpha_angle))/omega_ce) ;

chi=(z/vthe/( kx*cos(alpha_angle)  +  ky*sin(alpha_angle)    )  )  -  ( vde*sin(alpha_angle) / vthe );


return_residual= 1 - ...
    ((omega_pi*omega_pi/(kk*kk*vthi*vthi))*Zplasma_d(( z - ky*vdi)/(kk*vthi)  ))+...
    (2*omega_pe^2/(kk*kk*vthe*vthe) ) * (...
1+...
((z - ky*vde)/(( kx*cos(alpha_angle) +  ky*sin(alpha_angle)     )*vthe)*...
besselj(0,bb)*besseli(0,lam)*exp(-lam)*Zplasma(chi) ...
)...
);

if (isinf(return_residual)) return_residual=1000;end;
if (isnan(return_residual)) return_residual=1000;end;

output_args=[real(return_residual), imag(return_residual)];

end

