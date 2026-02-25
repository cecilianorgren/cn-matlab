function [ output_args ] = dispersion_relation_yoon2008(input_args,solve_parameterX,solve_parameterY,solve_parameterZ,physical_parameters)
% Attention: takes a complex argument as a (real(x), imag(x)) vector
% A special version, use Yoon/Lui (2008) parameters to study Buneman-LHDI
% interaction.

% Physical constants
kB = 1.381e-23; % J/K
e = 1.602e-19; % C

% Physical parameters
% input physical_parameters = [omega_pi,omega_pe,omega_ci,omega_ce,omega_lh,vthe,vthi,vdi,vde,R,Ln];
omega_pi = physical_parameters(1); % rad/s
omega_pe = physical_parameters(2); % rad/s
omega_ci = physical_parameters(3); % rad/s
omega_ce = physical_parameters(4); % rad/s
omega_lh = physical_parameters(5); % rad/s
vthe = physical_parameters(6); % m/s 
vthi = physical_parameters(7); % m/s 
vdi = physical_parameters(8);     % ion diamagnetic drift, m/s
vde = physical_parameters(9);    % electron diamagnetic drift, m/s
Rop = physical_parameters(10); % m
Ln = physical_parameters(11);                    % gradient length scale, m

% Guide field
alpha_angle = physical_parameters(12);           % need to be a small value (no guide field, see Yoon/Lui 2008 paper

if 0
disp(['ope=' num2str(omega_pe,'%.0f') ', opi=' num2str(omega_pi,'%.0f') ', oce=' num2str(omega_ce,'%.0f'),...
      ', oci=' num2str(omega_ci,'%.0f') ', olh=' num2str(omega_lh,'%.0f') ', vti=' num2str(vthi,'%.0f'),...
      ', vte=' num2str(vthe,'%.0f') ', vde=' num2str(vde,'%.0f') ', vdi=' num2str(vdi,'%.0f') ', rhop=' num2str(Rop,'%.0f'),...
      ', Ln=' num2str(Ln,'%.0f')])
end


% omega ?
z=input_args(1)+1i*input_args(2);

kx = solve_parameterX*omega_ce/vthe;  %omega_ci/vthi;%omega_ce/vthe;
ky = solve_parameterY*omega_ce/vthe;  %*omega_ci/vthi;%omega_ce/vthe;
kz = -1i*solve_parameterZ*omega_ce/vthe;%  - (1.2*solve_parameterZ*omega_ce/vthe); %*omega_ci/vthi;%*omega_ce/vthe;
kk = sqrt(kx.^2 + ky.^2 - kz.^2);

bb=((kx*sin(alpha_angle))-(ky*cos(alpha_angle)) )*vde*cos(alpha_angle)/omega_ce;
lam=((  vthe*vthe*  ( -kz.^2 + ( (kx*sin(alpha_angle)) - ( ky*cos(alpha_angle))   ).^2) )/(2*omega_ce*omega_ce)) + ((kz*vde*cos(alpha_angle))/omega_ce) ;
chi=(z/vthe/( kx*cos(alpha_angle)  +  ky*sin(alpha_angle)    )  )  -  ( vde*sin(alpha_angle) / vthe );

return_residual = 1 - ...
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

