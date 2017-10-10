function [ocs_kdir correlation corramp] = cn_tool_v2(B,E,Ne,flow,v,tries)
% Tries different propagation directions 
% and plots phi_B and phi_E.
% Gives back the different vectors x y z, 
% y is propagation direction,
% x is normal direction.

irf_units;
Ne=Ne*1e6;  % m^-3
row=size(E,1);
B0=mean(B(:,2:4),1); B0=sqrt(B0(1).^2+B0(2).^2+B0(3).^2); % B0 nT

% Resample
B=irf_resamp(B,E);

% Filter
fs=1/(E(2,1)-E(1,1));
b=irf_filt(B,flow,0,fs,5); b=b(:,2:4);
e=irf_filt(E,flow,0,fs,5); e=e(:,2:4);

% Set up coordinate systems
z=irf_norm(mean(B));
z=z(2:4);
x=irf_norm(cross(z,cross([1 1 0],z)));
y=cross(z,x);

% Get Phi_B    
facb=[b*x' b*y' b*z'];
phib=[B(:,1) repmat(facb(:,3),1,tries)*B0*1e-18/(Ne*Units.mu0*Units.e)];

theta=2*(0:pi/tries:pi-pi/tries)'; % half a turn, then take max(abs(corr))
fack=[cos(theta) sin(theta) 0*sin(theta)]; % evenly spread out in perp plane

% Propagation direction in original coordinate system
ocs_kdir=cn_m_trans(fack,[x;y;z],-1);

% Get e in "propagation direction", rows are time, col is direction
%ek=face*fack';
ek=e*ocs_kdir';
phie=irf_integrate([E(:,1) -detrend(ek*v,'constant')]); % [time phie1 phie2 ... phietries] 

%correlation=zeros(tries,2);
figure;
pp='plot';
if strcmp(pp,'plot')
if tries<20
    h=irf_plot(tries); 
    flag=1;
    for k=1:tries;
        irf_plot(h(k),{phie(:,[1 k+1]),phib},'comp');
        irf_legend(gca,{'phi_E','phi_B','amp'},[0.90,0.90])
    end
end
end
%return
correlation=zeros(tries,2);

for k=1:tries          
    % Check correlation
    corrmatrix=corrcoef(phib(:,2),phie(:,k+1));
    correlation(k,1)=corrmatrix(1,2);
    
    if 0
    if flag==1
        sumcorr=sum((phie(:,1+k)-phib(:,2)).^2);
        xcorr1=sumcorr0;
        flag=0;
    end
    end
    correlation(k,2)=sum((phie(:,1+k)-phib(:,2)).^2);
    % Check waveform
    %window=10;
    %[t0x(k),~,~,~,~]=cn_mycorr(PhiE(:,2),PhiB(:,2),window,450);
    %[t0x2(k) corrx2(:,k)]=cn_xcorr(PhiE,PhiB,window,'x');
    
   % irf_plot(h(k),{PhiE,PhiB},'comp')
   % irf_legend(h(k),{'phi_E','phi_B'},[0.90,0.90])
end
correlation(:,2)=correlation(:,2)/max(correlation(:,2));
1;
return
i1=find(correlation(:,2)==min(correlation(:,2)))
%Ek=[EAC(:,1) (x(i1,:)*EAC(:,(2:4))')'];
% Get Phi_E
%size(PhiE)
%    PhiE=irf_integrate(Ek);
%    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
%    PhiE=[PhiE(:,1) PhiE(:,2)-mean(PhiE(:,2))];
%h=irf_plot(2);

%irf_plot(h(1),{PhiE,PhiB},'comp')
irf_legend(h(1),{'phi_E','phi_B','amp'},[0.90,0.90])
%irf_zoom(h,'x',[PhiE(1,1) PhiE(end,1)])
i2=find(correlation(:,1)==max(correlation(:,1)))
Ek=[EAC(:,1) (x(i2,:)*EAC(:,(2:4))')'];
% Get Phi_E
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
    PhiE=[PhiE(:,1) PhiE(:,2)-mean(PhiE(:,2))];
%    irf_plot(h(2),{PhiE,PhiB},'comp')
%irf_legend(h(2),{'phi_E','phi_B','xcorr'},[0.90,0.90])

correlation;
corrx=correlation(:,1);
corramp=correlation(:,2);