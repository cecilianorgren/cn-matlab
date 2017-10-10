if 1 % the beam model
n  = [0.060 0.0025]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.500 0.060]; % keV
vd = [0 -4];
d  = [1 1]; % no loss cone
a1 = [1 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 0;
tbeam = toepoch([2007 08 31 10 17 38.2]);
end

in=2;
%%
h_ax=axes; hold(h_ax,'on')
[h,f,vtot]=whamp.plot_f(n(in),m(in),t(in),vd(in),d(in),a1(in),a2(in),pitchangles,plotoption);

for k=1:3
    density{k}=trapz(vtot(k,:),f(k,:))*2*pi; % /km^3   
end

disp([num2str(density{3}) ' particles/km^3 = ' num2str(density{3}*1e-15) ' cc'])
% v is m/s, f is s^3/m^6
%%
% testar att integrera f bara och ser om ja f?r tillbaka densiteten

% f?rst integrerar jag ?ver vinkelr?ta riktningen
    n1=trapz(vtot(2,:),f(2,:)); % s/m^4


%%
% integrate 

for k=1:3
for p=3:size(f,2)
    integral{k}(p)=trapz(vtot(k,1:p),f(k,1:p).*vtot(k,1:p)); % s/m^4    
end
    current{k}=trapz(vtot(k,:),f(k,:).*vtot(k,:)); % s/m^4
    density{k}=trapz(vtot(k,:),f(k,:)); % s/m^4   
end

for k=1:3; h(k) = subplot(3,1,k); end

subplot(3,1,1); [h,f,vtot]=whamp.plot_f(n(in),m(in),t(in),vd(in),d(in),a1(in),a2(in),pitchangles,plotoption);
plot(h(2),vtot(3,:),f(3,:)) 
plot(h(3),vtot(3,:),integral{3}) 

%%
[h,f,vtot]=whamp.plot_f(n(in),m(in),t(in),vd(in),d(in),a1(in),a2(in),pitchangles,2);
hold(gca,'on')
%plot(vtot(3,:),f(3,:))
plot(vtot,f*2)
hold(gca,'off')
%units=irf_units;
%diffe=(integral{3}-integral{1})
%je=diffe*units.e