function [correlation] = tool_velocity2(B0,Bz,PhiE,n,v,title_str)
% Finds best match in amplitude given, B0, dB_par, phi, propagation 
% direction implied, for specified n and v given as vectors.
% Returns a matrix of correlations and A, im and map to make a gif if
% visualization is toggled.
%
% phi = B0*dB_par/n_e*e*mu0
% phi = int(E)dt*v (dl=-vdt => -dl=vdt)
%

% Adjust title_str and add density and B0
title_str=[title_str,', ',num2str(B0,'%.f'),' nT'];

% Define constants
mu0=4*pi*1e-7;
n=n*1e6; % density in #/m^3
e=1.6e-19;

% Allocate correlations matrix 
nn=length(n);
nv=length(v);
correlation=zeros(nn,nv); % rows: n, cols: v

%ylims=[floor(min(PhiE(:,2))*max(v)/100) ceil(max(PhiE(:,2))*max(v)/100)]*100;

LHS=[PhiE(:,1) PhiE(:,2:end)*torow(v)]; % varying v
RHS=[Bz(:,1) Bz(:,2)*B0*1e-18/mu0/e*(1./torow(n))]; % varying n

for k=1:nn;
    for p=1:nv;
        correlation(k,p)=sum((LHS(:,1+p)-RHS(:,1+k)).^2);
    end
end

if 0
fig=figure(84);
set(gcf,'color','white'); % white background for figures (default is grey)
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
ind=0; n_ims=tries;
for k=1:tries
    ind=ind+1;
    %LHS=[PhiE(:,1) PhiE(:,2)*v(k)];
    %RHS=[Bz(:,1) Bz(:,2)*B0*1e-18/mu0/N/e];    
    correlation(k,1)=sum((LHS(:,1+k)-RHS(:,2)).^2);
    h=irf_plot({LHS(:,[1 1+k]),RHS},'comp'); 
    irf_legend(h,['v = ',num2str(v(k),'%0.f'),' km/s'],[0.02 0.95]);
    irf_zoom(h,'x',[LHS(1,1) LHS(end,1)]);
    set(h,'ylim',ylims); grid off;
    title(h,title_str)
    hold off;
    
    
     % collect frames    
    f=getframe(fig);
    A(:,ind)=f;
    if k==1, % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,n_ims) = 0;
    else
        im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
    end
end
end