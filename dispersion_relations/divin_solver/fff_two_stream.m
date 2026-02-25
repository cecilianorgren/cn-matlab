function [ output_args ] = fff_two_stream( input_args,solve_parameter)
% Attention: takes a complex argument as a (real(x), imag(x)) vector


% omegab1=10;
% omegab2=10;
% 
% vb0=0.2;
% vb1=0.4;
% vb2=0.8;
% %qq=solve_parameter*sqrt(mm);
% 
% kk=solve_parameter;
% z=input_args(1)+1i*input_args(2);
% 
% return_residual= 1 -( ((omegab1^2)/(z-kk*vb0)^2)/1836    ) - ((omegab1^2)/(z-kk*vb1)^2);%-((omegab2^2)/(z-kk*vb2)^2);




%% execute that in Matlab before running solver script.

% Te=1.5 KeV, Ti=2 KeV
vthe=sqrt(2)*1.62e7;            % electron thermal velocity 1, m/s
vthi=sqrt(2)*4.37e5;            % ion thermal velocity, m/s
vtheB=sqrt(2)*3.24e6;           % electron thermal velocity 2, m/s
vbeB=1.*vthe;                   % electron drift velocity 2, m/s

vtheBCOLD=3.24e4;
vbeCOLD=0;

vBe=-2.8*vthe                   % electron drift velocity 1, m/s


mproton=1.6726*(10^(-27)); % kg
qproton=1.6022*(10^(-19)); % C
melectron=9.10938291*(10^(-31)) %kg


omega_pe=1.38e4;                % electron plasma frequency 1, m/s
omega_pi=322;                   % ion plasma frequency, m/s
omega_peB=2820;                 % electron plasma frequency 2, m/s
omega_peBCOLD=0.5e4;

vSound=sqrt(((1*1*melectron*vthe*vthe)+(3*mproton*vthi*vthi))/mproton );

de= 299792458/omega_pe;
di= 299792458/omega_pi;

%length scale
lamD=vthi/omega_pi;
%%

z=input_args(1)+1i*input_args(2);

kx = solve_parameter/lamD;  %omega_ci/vthi;%omega_ce/vthe;


%kx = solve_parameter/di;  %omega_ci/vthi;%omega_ce/vthe;

return_residual=1-((omega_pi*omega_pi/(kx*kx*vthi*vthi )) *   (Zplasma_d( (z - kx*0 )/( kx*vthi )     )   )     ) ...
     - ((omega_pe*omega_pe/(kx*kx*vthe*vthe )) *   (Zplasma_d( (z - kx*vBe )/( kx*vthe )     )   )     );% ...
  %   - ((omega_peB*omega_peB/(kx*kx*vtheB*vtheB)) *   (Zplasma_d( (z - kx*vbeB )/( kx*vtheB )     )   )     ); %...
  %   - ((omega_peBCOLD*omega_peBCOLD/(kx*kx*vtheBCOLD*vtheBCOLD)) *   (Zplasma_d( (z - kx*vbeCOLD )/( kx*vtheBCOLD )     )   )     ) ...
  % ;

%return_residual= 1 - (omega_pi*omega_pi/(z*z)) - (omega_pe*omega_pe/(z*z)) - (omega_peB*omega_peB/((z-(kx*vbeB))^2 ));




if (isinf(return_residual)) return_residual=1000;end;
if (isnan(return_residual)) return_residual=1000;end;


% return_residual=(qq^2*epsilon) - ...
% Plasmaf_derivative((z-qq*uu)/qq) + ...
% (2/ttt)*( ...
% 1+ (( (1/sqrt(mm*ttt))*besselj(0,(qq*uu*ttt*cos(alpha_angle)*cos(alpha_angle)/mm)) *  (z+qq*uu*ttt)/(qq*sin(alpha_angle))    )*  ...
% besseli(0,(qq*qq*ttt*cos(alpha_angle)*cos(alpha_angle)/(2*mm)))*...
% exp((qq*qq*ttt*cos(alpha_angle)*cos(alpha_angle)/(2*mm)))*...
% Plasmaf( (1/sqrt(mm*ttt))*(z+qq*uu*ttt*sin(alpha_angle)*sin(alpha_angle))/(qq*sin(alpha_angle))   )...
% ) ...
% )



% Plasmaf_derivative((z-qq*uu)/qq)
% besselj(0,(qq*uu*ttt*cos(alpha_angle)*cos(alpha_angle)/mm)/mm)
% besseli(0,(qq*qq*ttt*cos(alpha_angle)*cos(alpha_angle)/(2*mm)))
% exp((qq*qq*ttt*cos(alpha_angle)*cos(alpha_angle)/(2*mm)))
% Plasmaf( (1/sqrt(mm*ttt))*(z+qq*uu*ttt*sin(alpha_angle)*sin(alpha_angle))/(qq*sin(alpha_angle))   )
% 
%  (1/sqrt(mm*ttt))*(z+qq*uu*ttt*sin(alpha_angle)*sin(alpha_angle))/(qq*sin(alpha_angle))  
%return_residual=Plasmaf(z)-solve_parameter;

output_args=[real(return_residual), imag(return_residual)];


end

