%% Beam
n   = [0.06e6 0.0025e6 0.0625e6 0 0 0 0 0 0 0];        % density
t   = [1.5 0.06 2.2 0 0 0 0 0 0 0];            % temperature [keV]
d   = [1 1 1 1 1 1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
a   = [1 1 1 1 1 1 1 1 1 1];           % T_perp/T_par
b   = [0.0 0.0 0 0 0 0 0 0 0 0];        % loss cone parameter
mass= [0 0 1 0 0 0 0 0 0 0];        % mass
vd  = [0 -4 0 0 0 0 0 0 0 0];         % V_drift/V_term
fce = 0.670; % kHz
pzl = 1;  
pitchangle = [0 45 90 135 180];
plotoption = 1; % f vs E
excl = []; 
n(excl)=0;
%[h,f,Etot] = whamp.plot_f(n,mass,t,vd,d,a,b,pitchangle,plotoption);

%% No beam
n   = [0.0625e6 0 0.0625e6 0 0 0 0 0 0 0];        % density
t   = [1.5 0.060 2.2 0 0 0 0 0 0 0];            % temperature [keV]
d   = [1 1 1 1 1 1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
a   = [1 1 1 1 1 1 1 1 1 1];           % T_perp/T_par
b   = [0 0 0 0 0 0 0 0 0 0];        % loss cone parameter
mass= [0 0 1 0 0 0 0 0 0 0];        % mass
vd  = [0 0 0 0 0 0 0 0 0 0];         % V_drift/V_term
fce = 0.670; % kHz
pzl = 1;  
pitchangle = [0 45 90 135 180];
plotoption = 1; % f vs E
excl = [2]; 
n(excl)=0;
%%
if 0
%[h,f,Etot] = whamp.plot_f(n,mass,t,vd,d,a,b,pitchangle,plotoption);
n   = [0.0625e6 0 0.0625e6 0 0 0 0 0 0 0];        % density
t   = [2.2 0.060 1.5 0 0 0 0 0 0 0];            % temperature [keV]
d   = [1 1 1 0 0 0 0 0 0 0];            % loss cone parameter, 1.0 = no loss cone)
a   = [1 1 1 0 0 0 0 0 0 0];           % T_perp/T_par
b   = [0 0 0 0 0 0 0 0 0 0];        % loss cone parameter
mass= [1 0 0 0 0 0 0 0 0 0];        % mass
vd  = [0 0 0 0 0 0 0 0 0 0];         % V_drift/V_term
fce = 0.670; % kHz
pzl = 1;  
pitchangle = [0 45 90 135 180];
plotoption = 1; % f vs E
inc = [1 3]; 
%[h,f,Etot] = whamp.plot_f(n,mass,t,vd,d,a,b,pitchangle,plotoption);
end
%%
myformat= [repmat('%.6f ',1,10) '\n'];

% open the file with permission to append
txtModel = '20130831-beam.txt';
fid = fopen(['/Users/Cecilia/Courses/PlasmaWaves/Whamp/Models/' txtModel],'w');

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

%%
% now run whamp and save data
%%
cd      = '/Users/Cecilia/irfu-matlab/+whamp/';
txtFile = '20130831-alfven.txt';
txtFile = '20130831-test.txt';
fname1  = ['/Users/Cecilia/Courses/PlasmaWaves/Whamp/Results/' txtFile];
file1   = load(fname1);

fig=figure(1);
set(fig,'defaultAxesFontUnits','pixels');
set(fig,'defaultAxesFontSize',14);
set(fig,'defaultTextFontSize',14);

[p,z,f,fim]=m2xyz(file1); 
their.p=p;their.z=z;their.f=f;their.fim=fim;
inst=[];%find(fim>0);
surf(log10(p(:)),log10(z(:)),10.^f,10.^fim) ; hold on;

title([txtModel ' -> ' txtFile])
xlabel(gca,'log_{10}(p)','fontsize',14);
ylabel(gca,'log_{10}(z)','fontsize',14);
zlabel(gca,'f/f_{cp}','fontsize',14);
cax=colorbar;
%ylabel(cax,'f_i/f_{pp}')
set(gcf,'paperpositionmode','auto');