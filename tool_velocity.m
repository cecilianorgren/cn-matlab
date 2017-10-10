function [correlation,vis] = tool_velocity(B0,Bz,PhiE,n,v,title_str,visualize)
% IRF_MATCH_PHIB_V
%
% Used together with irf_match_phib_dir.m. Finds best match in amplitude 
% given, B0, dB_par, phi, propagation direction implied, for specified n 
% and v given as vectors. Returns a matrix of correlations and A, im and 
% map to make a gif if visualization is toggled.
%
% phi = B0*dB_par/n_e*e*mu0
% phi = int(E)dt*v (dl=-vdt => -dl=vdt)
%

% Define constants
mu0=4*pi*1e-7;
n=n*1e6; % density in #/m^3
e=1.6e-19;

% Allocate correlations matrix 
nn=length(n);
nv=length(v);
correlation=zeros(nn,nv); % rows: n, cols: v

% Setup potentials 
LHS=[PhiE(:,1) PhiE(:,2:end)*torow(v)]; % depends on v
RHS=[Bz(:,1) Bz(:,2)*B0*1e-18/mu0/e*(1./torow(n))]; % depends on n

% Get correlation
for k=1:nn;
    for p=1:nv;
        correlation(k,p)=sum((LHS(:,1+p)-RHS(:,1+k)).^2);
    end
end

vis=[];
if visualize 
    switch nn 
        case 1 % make a gif with different velocities
            % Adjust title_str and add density and B0
            if ~exist('title_str','var') 
                title_str=[];
            end
            title_str=[title_str,', ',num2str(B0,'%.f'),' nT, ',num2str(n,'%.2f'), 'cc'];
            ylims=[floor(min(PhiE(:,2))*max(v)/100) ceil(max(PhiE(:,2))*max(v)/100)]*100;

            fig=figure;
            set(gcf,'color','white'); % white background for figures (default is grey)
            set(gcf,'defaultAxesFontSize',14);
            set(gcf,'defaultTextFontSize',14);
            ind=0;
            for k=1:nv
                ind=ind+1;
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
                    im(1,1,1,nv) = 0;
                else
                    im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
                end
            end
            vis.A=A;
            vis.im=im;
            vis.map=map;
        otherwise
            pcolor(v,n,log10(correlation));
            vis.h=gca;
            vis.hc=colorbar;
    end
else
    vis.A=[]; vis.im=[]; vis.map=[];
end