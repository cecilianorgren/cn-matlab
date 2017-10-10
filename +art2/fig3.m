% This script produces an overview plot of the spacecraft configuration for
% the time of the observation. It uses the same 2D coordinate system as is
% used to present the electric field. NO! That's not the case. The electric
% field y-component is the one that lies in the spin plane. 
% Maybe there's no problem since they are probably very close.

cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields

if ~exist('ind','var'); ind = 5; end
tint = tint{highQuality(ind)};

tint=tt;

%d00 = irf_tlim(d0,tint);
x=cn.mean(irf_tlim(dp1,tint),1); % SP
y=cn.mean(irf_tlim(dp2,tint),1); % BxSP
z=cn.mean(irf_tlim(dp3,tint),1); % B

x3 = x*(0.5);
y3 = y*(0.5);
z3 = z*(0.5);
x4 = x*(-0.5);
y4 = y*(-0.5);
z4 = z*(-0.5);

%figure(23)
set(gcf,'defaultAxesFontSize',16);
set(gcf,'defaultTextFontSize',16);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

h = plot3(x3,y3,z3,'go',x4,y4,z4,'bo');
%h = plot3(x3,y3,z3,'go',x4,y4,z4,'bo');
set(h,'markersize',14,'linewidth',2,'marker','square')%,'ydata',y*1.2*[-1 1],'zdata',z*1.2*[-1 1])

axis equal;
xlabel('\perp 2,sp')
ylabel('\perp 1')
zlabel('||,B')
view([-1 0 0])

text(x3,y3,z3,'      C3')
text(x4,y4,z4,'      C4')
axis([x*0.8*[-1 1] y*0.8*[-1 1] z*0.7*[-1 1]])



if 1 % 2D
h = plot(y3,z3,'go',y4,z4,'bo');
set(h,'markersize',14,'linewidth',2,'marker','square')%,'ydata',y*1.2*[-1 1],'zdata',z*1.2*[-1 1])

axis equal;
title(['\Delta \perp 2,sp = ' num2str(x,'%0.0f') ' km'])
xlabel('\perp 1')
ylabel('||,B')
view([0 0 -1])

text(y3-6,z3,'C3      ')
text(y4,z4,'     C4')
%axis([y*0.87*[-1 1] z*0.7*[-1 1]]) 
axis([[-15 15] [-19 19]]) 
set(gca,'xtick',[ -10 -5 0 5 10],'ytick',[-15 -10 -5 0 5 10 15])
end


%p = cn.mean(irf_tlim(dp1,tint),1);
%dm = irf_tappl(poss,'*0.5');

%c_eval('rC?(1) = irf_dot(dm,first); rC?(2) = irf_dot(dm,second); rC?(3) = irf_dot(dm,third);',3:4);
%rC4 = -dm;

% rC34 = [rC3;rC4];
% xC34 = rC34(:,1);
% yC34 = rC34(:,2);
% zC34 = rC34(:,3);
% 
% plot3(rC3(1),rC3(1),rC3(1),'go'); hold on;
% plot3(rC4(1),rC4(1),rC4(1),'bo'); hold off;
%d1 = irf_dot(poss,irf_tlim(first,tint));
%d2 = irf_dot(poss,irf_tlim(second,tint));
%d3 = irf_dot(poss,irf_tlim(third,tint));

