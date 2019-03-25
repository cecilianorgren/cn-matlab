% Illustrate fourier transform

syms t v w x f(x)
v = diff(x,t);
f =  exp(-(x-v*t).^2);

fourier_transform = fourier(f,x,t)%   returns   pi^(1/2)*exp(-t^2/4)

phi = @(t,x,v) exp(-(x-v*t).^2);

%% Single pulse
syms t w x phi E
v = 1;
phi =  exp(-(x-v*t)^2);
E = diff(phi,x);
ft = fourier(E,w,x);
mf_ft = matlabFunction(ft);
mf_phi = matlabFunction(phi);
mf_E = matlabFunction(E);


x_vec = linspace(-4,4,100);
plot(x_vec,mf_phi(x_vec),...
     x_vec,mf_E(x_vec),...
     x_vec,imag(mf_ft(x_vec)))
   
%% Fast fourier transform of travelling solitary pulse
v = 1;
x = 0;
syms t x phi
phi =  exp(-(x-v*t)^2);
E = diff(phi,x);

mf_phi = matlabFunction(phi);
mf_E = matlabFunction(E);

Fs = 1000;           % Sampling frequency
t = -0.5:1/Fs:0.5;  % Time vector 
L = length(t);      % Signal length

n = 2^nextpow2(L);
f = Fs*(0:(n/2))/n;
phi_of_t = mf_phi(t,0);
fft_phi = fft(phi_of_t,n);
P_phi = abs(fft_phi/n);

E_of_t = mf_E(t,0);
fft_E = fft(E_of_t,n);
P_E = abs(fft_E/n);

mf_base = @(t,f) cos(f.*t);
[T,F] = meshgrid(t,f(1:n/2+1));
[T,P] = meshgrid(t,P_E(1:n/2+1));
mf_BASE = mf_base(T,F);

%mf_BASE_P  = mf_BASE.*P_E;
%mf_BASE_F_cumsum = cumsum(mf_BASE_P,1);

% Plot
nrows = 4;
ncols = 2;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

isub = 1;
hca = h(isub); isub = isub + 1;
plot(hca,t_vec,mf_E(t_vec,0))

hca = h(isub); isub = isub + 1;
plot(hca,f,P_E(1:n/2+1))

hca = h(isub); isub = isub + 1;
plot(hca,f,real(fft_E(1:n/2+1)))

hca = h(isub); isub = isub + 1;
plot(hca,f,imag(fft_E(1:n/2+1)))


hca = h(isub); isub = isub + 1;
imagesc(hca,t,f,(mf_BASE))

hca = h(isub); isub = isub + 1;
imagesc(hca,t,f,cumsum(mf_BASE,1))

hca = h(isub); isub = isub + 1;
imagesc(hca,t,f,cumsum(mf_BASE.*F,1))


%%

for ifreq = 1:2:size(mf_BASE_F,1)
  plot(t,mf_BASE_F_cumsum(ifreq,:));
  pause(0.1)
end
