% run whamp for electron hole article #2

% values from figure 8distsandefield_1-8.png (expect ions)
% hot electrons/cold drifting electrons/hot ions
n =  [0.055e6 0.0075e6 (0.055+0.0075)*1e6 0 0 0 0 0 0 0];  % density
t =  [1.6 0.060 2 0 0 0 0 0 0 0];  % temperature [keV]
d =  [1 1 1 1 1 1 1 1 1 1];        % loss cone parameter, 1.0 = no loss cone)
a =  [1 1 1 1 1 1 1 1 1 1];        % T_perp/T_par
b =  [0.0 0.0 0 0 0 0 0 0 0 0];    % loss cone parameter
mass=[0 0 1 0 0 0 0 0 0 0];        % mass
vd = [0 -3.5 0 0 0 0 0 0 0 0];     % V_drift/V_term
vd = [0 -7.5 0 0 0 0 0 0 0 0];     % V_drift/V_term
fce = 0.7;                         % electron gyrofrequency in kHz
pzl=0;                             % /1 -log scale, 0 - linear scale /   
%opens tr
%pa=[0 90 180];        % pitch angeles
%plotoption=[1];
%titleoption=[0];

if 0
%%

%opens tr
%pa=[0 90 180];        % pitch angeles
%plotoption=[1];
%titleoption=[0];
units=irf_units;
vth=sqrt(2*units.e*t*1000./mass);
ffcp=[1 3];
z=[1e-3 1];
p=[1e-4 1e-1];
end
%%
myformat= [repmat('%.6f ',1,10) '\n'];

% open the file with permission to append
filePath = '/Users/Cecilia/Whamp/Models/';
fileName = 'mod_buneman.txt';
fid = fopen([filePath fileName],'w');

% write values at end of file
fprintf(fid, myformat,n);
fprintf(fid, myformat,t);
fprintf(fid, myformat,d);
fprintf(fid, myformat,a);
fprintf(fid, myformat,b);
fprintf(fid, myformat,mass);
fprintf(fid, myformat,vd);
fprintf(fid, '%.3f \n',fce);
fprintf(fid, '%1.0f \n',pzl);
fprintf(fid, myformat,vd);

% close the file 
fclose(fid);
eval(['type ' filePath fileName])

%% what to enter in whamp
% p - perpendicular wavenumber
% p = 0 - parallel propagation
% we want to find surface at 0.01 omega_pe
fpe = 9*sqrt(n(3)*1e-6);
units = irf_units;
vth = sqrt(2*units.e*t(1)*1e3/units.me);
f_real=0.01*fpe;
f = f_real/fpe; % in units of fpe 
v_real = 1000; % km/s
l_par = v_real/f_real;
z = vth/l_par/fce; % ~700

% so we could take
f = [0.001 0.1];
z = [10 2000];
%% run whamp
%% visualize
cd = '/Users/Cecilia/irfu-matlab/+whamp/';
fname = [filePath fileName];
file = load(fname);
%close 1
fig=figure(1);
set(fig,'defaultAxesFontUnits','pixels');
set(fig,'defaultAxesFontSize',14);
set(fig,'defaultTextFontSize',14);
%test=/Users/Cecilia/Courses/PlasmaWaves/Whamp/Results/res1.txt;
% subplot(2,2,1)
% [p,z,f,fim]=m2xyz(file1); 
% surf(log10(p),log10(z),f,fim) 
% colorbar()
% xlabel('p');
% ylabel('z');
% zlabel('f');
% subplot(2,2,2)
% [p,z,f,fim]=m2xyz(file2); 
% surf(log10(p),log10(z),f,fim) 
% colorbar()
% xlabel('p');
% ylabel('z');
% zlabel('f');

%subplot(2,2,3)
[p,z,f,fim]=m2xyz(file); 
their.p=p;their.z=z;their.f=f;their.fim=fim;
inst=[];%find(fim>0);
surf(log10(p(:)),log10(z(:)),f,fim) ; hold on;

%[p,z,f,fim]=m2xyz(file2); 
%surf(log10(p),log10(z),f,fim) ; hold on;
%[p,z,f,fim]=m2xyz(file3); 
%surf(log10(p),log10(z),f,fim) ; hold on;
%[p,z,f,fim]=m2xyz(file4); 
%surf(log10(p),log10(z),f,fim) ; hold off;
xlabel(gca,'log_{10}(p)','fontsize',14);
ylabel(gca,'log_{10}(z)','fontsize',14);
zlabel(gca,'f/f_{pp}','fontsize',14);
cax=colorbar;
ylabel(cax,'f_i/f_{pp}')
set(gcf,'paperpositionmode','auto');
%% try loading electric fields

%fn=fopen('/Users/Cecilia/Courses/PlasmaWaves/Whamp/Results/iacw_2e.txt');
%fl=textscan(fn,'%*s%f%f%*s%f%f%*s%f%f\n')
%
fn=fopen('/Users/Cecilia/Courses/PlasmaWaves/Whamp/Results/whistler_pzfe3.txt');
text_format='%*s%f%f%*s%f%f%*s%f%f';
text_format2='%s';

data=textscan(fn,text_format2,'delimiter','\n');
%isEmtpyLine= cellfun(@(line)isempty(line),regexp(data{1},'^[- ]*$'))'
%data=data(isEmptyLine);
%%
nii=length(data{1});
empty=zeros(1,nii);
p0=zeros(nii,1);
z0=zeros(nii,1);
f0=zeros(nii,1);
fim0=zeros(nii,1);
Ex0=zeros(nii,2);
Ey0=zeros(nii,2);
Ez0=zeros(nii,2);
for ii=1:nii
    if ~strcmp(data{1}{ii},'') 
        empty(ii)=1;
        p0(ii)=str2double(data{1}{ii}(1:9));
        z0(ii)=str2double(data{1}{ii}(16:24));
        f0(ii)=str2double(data{1}{ii}(29:41));
        fim0(ii)=str2double(data{1}{ii}(43:51));
        Ex0(ii,1:2)=[str2num(data{1}{ii}(58:64)) str2num(data{1}{ii}(66:72))];
        Ey0(ii,1:2)=[str2num(data{1}{ii}(78:84)) str2num(data{1}{ii}(86:92))];
        Ez0(ii,1:2)=[str2num(data{1}{ii}(98:104)) str2num(data{1}{ii}(106:112))];
    end
end
empty=logical(empty);
p0=p0(empty);z0=z0(empty);f0=f0(empty);fim0=fim0(empty);
Ex0=Ex0(empty,:);Ey0=Ey0(empty,:);Ez0=Ez0(empty,:);

fille=[p0 z0 f0 fim0];
[p,z,f,fim]=m2xyz([p0 z0 f0 fim0]);
[~,~,Ex,Exi]=m2xyz([p0 z0 Ex0(:,1) Ex0(:,2)]);
[~,~,Ey,Eyi]=m2xyz([p0 z0 Ey0(:,1) Ey0(:,2)]);
[~,~,Ez,Ezi]=m2xyz([p0 z0 Ez0(:,1) Ez0(:,2)]);
own.p=p;own.z=z;own.f=f;own.fim=fim;
surf(log10(p),log10(z),f,Ez) ; hold on;
colorbar
%%
%close
Exx=Ex+1i*Exi;
Eyy=Ey+1i*Eyi;
Eperp=sqrt(Exx.*conj(Exx)+Eyy.*conj(Eyy));
Epolar=-2*imag(Exx.*conj(Eyy))./(Eperp.*Eperp);
surf(log10(p),log10(z),f,Ex./Ez) ; hold on;
colorbar
%%
%[p,z,f,Ex]=m2xyz(file2); 
%[p,z,f,Ex]=m2xyz(file2); 
surf(log10(p),log10(z),f,Epolar) ; hold on;
%[p,z,f,fim]=m2xyz(file3); 
%surf(log10(p),log10(z),f,fim) ; hold on;
%[p,z,f,fim]=m2xyz(file4); 
%surf(log10(p),log10(z),f,fim) ; hold off;
xlabel(gca,'log_{10}(p)','fontsize',14);
ylabel(gca,'log_{10}(z)','fontsize',14);
zlabel(gca,'f/f_{ce}','fontsize',14);
cax=colorbar;
ylabel(cax,'Polarization','fontsize',14)
set(gcf,'paperpositionmode','auto');
caxis(gca,[-1 1])