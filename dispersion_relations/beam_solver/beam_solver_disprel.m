function [ output_args ] = beam_solver_disprel(input_args,solve_parameter,fv)
% Attention: takes a complex argument as a (real(x), imag(x)) vector
%% execute that in Matlab before running solver script.
% fv = [omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2];
omega_pi = fv(1);
omega_pe1=fv(2);
omega_pe2=fv(3);
vthi=fv(4);
vthe1=fv(5);
vthe2=fv(6);
vdi=fv(7);
vde1=fv(8);
vde2=fv(9);
lamD=fv(10);
% Common parameters
%B = 25;
%n = 0.04;
%R = 0.20;
%n1= n*(1-R);
%n2= n*R;
%n1 = 0.04;
%n2 = n1*0.05;
%n = n1+n2;
%no = 0;
%Te1 = 1600; TeK1 = Te1*11604; % eV -> K
%Te2 = 60;   TeK2 = Te2*11604; % eV -> K
%Ti = 2000;  TiK  = Ti*11604; % eV -> K

%lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m

% Ions
%omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
%vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 

% Background electrons
% omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
% vthe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Vte'); % m/s 
% Ld1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Ld'); % m/s 
% vde1 = 0; % m/s
% 
% % Beam electrons
% omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
% vthe2 = irf_plasma_calc(B,n,no,Te2,Ti,'Vte'); % m/s 
% vde2 = S*vthe1;
% 
% % Other parameters
% omega_bune = omega_pe1^(1/3)*omega_pi^(2/3)*16^(-1/3);
% lamD = (1*vthe1/omega_pe1+1*vthe2/omega_pe2)/sqrt(2);

%%

z=input_args(1)+1i*input_args(2);

kx = solve_parameter/lamD; 

% Dispersion relation
return_residual = 1 - ((omega_pi*omega_pi/(kx*kx*vthi*vthi )) *   (Zplasma_d( (z - kx*vdi )/( kx*vthi )     )   )     ) ...
                    - ((omega_pe1*omega_pe1/(kx*kx*vthe1*vthe1 )) *   (Zplasma_d( (z - kx*vde1 )/( kx*vthe1 )     )   )     ) ...
                    - ((omega_pe2*omega_pe2/(kx*kx*vthe2*vthe2 )) *   (Zplasma_d( (z - kx*vde2 )/( kx*vthe2 )     )   )     );


if (isinf(return_residual)) return_residual=1000;end;
if (isnan(return_residual)) return_residual=1000;end;
output_args=[real(return_residual), imag(return_residual)];

end

function out_f=Zplasma_d(in_f)
%out_f=i*sqrt(pi)*exp(-in_f.*in_f).*(-2*in_f).*erfc_complex(-i*in_f)-2;
out_f=-2*(1+(in_f.*Zplasma(in_f)));
end

function out_f=Zplasma(in_f)
out_f =  1i*sqrt(pi)*Faddeeva(in_f,64);
end

function w = Faddeeva(z,N)
% FADDEEVA   Faddeeva function
%   W = FADDEEVA(Z) is the Faddeeva function, aka the plasma dispersion
%   function, for each element of Z. The Faddeeva function is defined as:
%
%     w(z) = exp(-z^2) * erfc(-j*z)
%
%   where erfc(x) is the complex complementary error function.
%
%   W = FADDEEVA(Z,N) can be used to explicitly specify the number of terms
%   to truncate the expansion (see (13) in [1]). N = 16 is used as default.
%
%   Example:
%       x = linspace(-10,10,1001); [X,Y] = meshgrid(x,x); 
%       W = faddeeva(complex(X,Y)); 
%       figure; 
%       subplot(121); imagesc(x,x,real(W)); axis xy square; caxis([-1 1]); 
%       title('re(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%       subplot(122); imagesc(x,x,imag(W)); axis xy square; caxis([-1 1]);
%       title('im(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%
%   Reference:
%   [1] J.A.C. Weideman, "Computation of the Complex Error Function," SIAM
%       J. Numerical Analysis, pp. 1497-1518, No. 5, Vol. 31, Oct., 1994 
%       Available Online: http://www.jstor.org/stable/2158232

if nargin<2, N = []; end
if isempty(N), N = 32; end

w = zeros(size(z)); % initialize output

%%%%%
% for purely imaginary-valued inputs, use erf as is if z is real
idx = real(z)==0; %
w(idx) = exp(-z(idx).^2).*erfc(imag(z(idx)));

if all(idx), return; end
idx = ~idx;

%%%%%
% for complex-valued inputs

% make sure all points are in the upper half-plane (positive imag. values)
idx1 = idx & imag(z)<0;
z(idx1) = conj(z(idx1));

M = 2*N;
M2 = 2*M;
k = (-M+1:1:M-1)'; % M2 = no. of sampling points.
L = sqrt(N/sqrt(2)); % Optimal choice of L.

theta = k*pi/M;
t = L*tan(theta/2); % Variables theta and t.
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f]; % Function to be transformed.
a = real(fft(fftshift(f)))/M2; % Coefficients of transform.
a = flipud(a(2:N+1)); % Reorder coefficients.

Z = (L+1i*z(idx))./(L-1i*z(idx));
p = polyval(a,Z); % Polynomial evaluation.
w(idx) = 2*p./(L-1i*z(idx)).^2 + (1/sqrt(pi))./(L-1i*z(idx)); % Evaluate w(z).

% convert the upper half-plane results to the lower half-plane if necesary
w(idx1) = conj(2*exp(-z(idx1).^2) - w(idx1));
end