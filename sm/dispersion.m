% this is for Octave users
set (0, 'defaultaxesposition', [0.15, 0.15, 0.795, 0.795])
%
% beam velocity
v0 = 0.2
% note that background density is 1 and each beam has 1/2 denisty
% plasma frequency
wp = 1.0/sqrt(2)
% wave number
k = 0:.01:5.0;

% solution of the dispersion relation
omega = -i*(v0^2*k.^2 + wp^2 - wp*(4*k.^2*v0.^2 + wp^2).^(1/2)).^(1/2);
plot(k,omega,'LineWidth',3)
xlabel('k')
ylabel('omega')
grid on

omega_most_unstable = max(omega)
k_most_unstable = k(find(omega==omega_most_unstable))
