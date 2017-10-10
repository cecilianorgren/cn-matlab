%% load Daniels data
load('/Users/Cecilia/Research/EH2/Daniel/wavespeeds.mat')
yES = wavespeeds.vESwaves./wavespeeds.vtrESwaves;
xES = wavespeeds.vESwaves./wavespeeds.vthESwaves;
yEH = wavespeeds.veholes./wavespeeds.vtreholes;
xEH = wavespeeds.veholes./wavespeeds.vtheholes;

%% Make figure with differetn combinations of vph, vT and vtebg

plot3(wavespeeds.veholes,wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',...
      wavespeeds.vESwaves,wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx')
xlabel('v_{ph}')
ylabel('v_{te,bg}')
zlabel('v_{v_T}')
%%
nPl = 6;
for kk = 1:nPl; h(kk) = subplot(2,3,kk); end
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes,wavespeeds.vtheholes,'ko',wavespeeds.vESwaves,wavespeeds.vthESwaves,'rx')
xlabel(hca,'v_{ph}'); ylabel(hca,'v_{te,bg}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx')
xlabel(hca,'v_{te,bg}'); ylabel(hca,'v_{T}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtreholes,wavespeeds.veholes,'ko',wavespeeds.vtrESwaves,wavespeeds.vESwaves,'rx')
xlabel(hca,'v_{T}'); ylabel(hca,'v_{ph}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtreholes./wavespeeds.vtheholes,wavespeeds.veholes,'ko',wavespeeds.vtrESwaves./wavespeeds.vthESwaves,wavespeeds.vESwaves,'rx')
xlabel(hca,'v_{T}/v_{te,bg}'); ylabel(hca,'v_{ph}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes./wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',wavespeeds.vESwaves./wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx')
xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{T}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes./wavespeeds.vtheholes,wavespeeds.veholes./wavespeeds.vtreholes,'ko',...
         wavespeeds.vESwaves./wavespeeds.vthESwaves,wavespeeds.vESwaves./wavespeeds.vtrESwaves,'rx')
xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{ph}/v_{T}')

for kk = 1:6;
    set(h(kk),'yscale','log','xscale','log')
end
%%
plot3(wavespeeds.veholes,wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',...
      wavespeeds.vESwaves,wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx')
xlabel('v_{ph}')
ylabel('v_{te,bg}')
zlabel('v_{v_T}')
%% Plot 'theoretical' data + the observations
nPl = 8;
for kk = 1:nPl; h(kk) = subplot(3,3,kk); end
isub = 1;

plys = 1:floor(size(vph,2)/3); plxs = 1:size(vph,1); % small Rs
plym = ceil(size(vph,2)/3):floor(size(vph,2)*2/3); plxm = 1:size(vph,1); % middle Rs
plyh = ceil(size(vph,2)*2/3):size(vph,2); plxh = 1:size(vph,1); % big Rs

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes,wavespeeds.vtheholes,'ko',...
         wavespeeds.vESwaves,wavespeeds.vthESwaves,'rx',...
         vph,ones(size(vph))*vthebg,'gs')
xlabel(hca,'v_{ph}'); ylabel(hca,'v_{te,bg}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',...
         wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx',...
         ones(size(vph))*vthebg,vT,'gs',...
         ones(size(vph))*vthebg,vphi,'bs')
xlabel(hca,'v_{te,bg}'); ylabel(hca,'v_{T}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtreholes,wavespeeds.veholes,'ko',...
         wavespeeds.vtrESwaves,wavespeeds.vESwaves,'rx',...
         vT(plxs,plys),vph(plxs,plys),'cs',...
         vphi(plxs,plys),vph(plxs,plys),'cp',...
         vT(plxm,plym),vph(plxm,plym),'ms',...
         vphi(plxm,plym),vph(plxm,plym),'mp',...
         vT(plxh,plyh),vph(plxh,plyh),'bs',...
         vphi(plxh,plyh),vph(plxh,plyh),'bp')
xlabel(hca,'v_{T}'); ylabel(hca,'v_{ph}')

if 0
    %%
    h=subplot(1,1,1);
    fa = 3.5;
    isub = 1;
    hca = h(isub); isub = isub + 1;
    plot(hca,vT(plxs,plys),vph(plxs,plys),'cs',...
             vphi(plxs,plys)/fa,vph(plxs,plys),'rp',...
             vT(plxm,plym),vph(plxm,plym),'ms',...
             vphi(plxm,plym)/fa,vph(plxm,plym),'gp',...
             vT(plxh,plyh),vph(plxh,plyh),'bs',...
             vphi(plxh,plyh)/fa,vph(plxh,plyh),'yp',...
             wavespeeds.vtreholes,wavespeeds.veholes,'ko',...
             wavespeeds.vtrESwaves,wavespeeds.vESwaves,'kx')
    xlabel(hca,'v_{T}'); ylabel(hca,'v_{ph}')
    set(h(1),'yscale','lin','xscale','log')
end

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtreholes./wavespeeds.vtheholes,wavespeeds.veholes,'ko',...
         wavespeeds.vtrESwaves./wavespeeds.vthESwaves,wavespeeds.vESwaves,'rx',...
         vT/vthebg,vph,'gs',...
         vphi/vthebg,vph,'bs')
xlabel(hca,'v_{T}/v_{te,bg}'); ylabel(hca,'v_{ph}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes./wavespeeds.vtheholes,wavespeeds.vtreholes,'ko',...
         wavespeeds.vESwaves./wavespeeds.vthESwaves,wavespeeds.vtrESwaves,'rx',...
         vph/vthebg,vT,'gs',...
         vph/vthebg,vphi,'bs')
xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{T}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.veholes./wavespeeds.vtheholes,wavespeeds.veholes./wavespeeds.vtreholes,'ko',...
         wavespeeds.vESwaves./wavespeeds.vthESwaves,wavespeeds.vESwaves./wavespeeds.vtrESwaves,'rx',...
         vph/vthebg,vph./vT,'gs',...
         vph/vthebg,vph./vphi,'bs')
xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{ph}/v_{T}')

hca = h(isub); isub = isub + 1;
plot(hca,wavespeeds.vtheholes,wavespeeds.veholes./wavespeeds.vtreholes,'ko',...
         wavespeeds.vthESwaves,wavespeeds.vESwaves./wavespeeds.vtrESwaves,'rx',...
         ones(size(vph))*vthebg,vph./vT,'gs',...
         ones(size(vph))*vthebg,vph./vphi,'bs')
xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{T}')

for kk = 1:6;
    set(h(kk),'yscale','log','xscale','log')
    set(h(kk),'yscale','lin','xscale','lin')
end