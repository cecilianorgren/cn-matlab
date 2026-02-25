% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

% Distribution properties 
gridding='test_less';
switch gridding    
    case 'kS' 
        Te1 = 300;
        Te2 = Te1/25;
        Ti = Te1;
        S = 0.3:0.05:0.8;
        R = [0.2];  
        n = 1;
    case 'test' 
        Te1 = 300;
        Te2 = Te1/25;
        Ti = Te1;
        S = 0.3:0.1:0.7;
        R = [0.05:0.2:0.95];  
        n = 1;
    case 'test_less' 
        Te1 = 5000;
        Te2 = Te1/100;
        Ti = Te1/11;
        S = 0.3:0.2:0.7; S = [0.35];
        R = [0.1:0.3:0.9]; R = [0.99];
        n = 5;
    case 'full'
        Te1 = 300;
        Te2 = Te1/25;
        Ti = Te1;
        S = 0.1:0.05:1.5;
        R = [0.01 0.05:0.05:0.95 0.99];        
        n = 1;
end

n = n*1e6;
ne1 = n*(1-R);
ne2 = n*R;
ni = n;

% Physical variables
vith = sqrt(2*qe*Ti/mi);
veth1 = sqrt(2*qe*Te1/me);
veth2 = sqrt(2*qe*Te2/me);
vde1 = veth1*0.0;
vde2 = veth1*S;

wpe = sqrt((n)*qe^2/(me*eps));
wpe1 = sqrt(ne1*qe^2/(me*eps));
wpe2 = sqrt(ne2*qe^2/(me*eps));
wpi = sqrt(ni*qe^2/(mi*eps));
Ld = sqrt(veth2.^2)./wpe/sqrt(2);

% wavenumbers
nk = 100;
kvec = linspace(1e-5,5,nk)/Ld;
dk2 = kvec(2)-kvec(1);

nR = numel(R);
nv = numel(S);
nk = numel(kvec);
nTe1 = numel(Te1);
nTe2 = numel(Te2);
nTi = numel(Ti);

% Initialize variables
wr = cell(nR,nv);
wi = cell(nR,nv);
vph = cell(nR,nv);
residual = cell(nR,nv);
wimax = nan(nR,nv);
wrmax = nan(nR,nv);
kmax = nan(nR,nv);
vphmax = nan(nR,nv);


%% Loop starts here
hca = subplot(1,1,1); hold(hca,'on')
tic
for doposwi = [1];
for iR = 1:nR
    for iS = 1:nv
        disp(['R=' num2str(R(iR)) ' S=' num2str(S(iS))])
        disp(sprintf('wpe = %s, vde = %s, veth = %s, wpi = %s, vith = %s',wpe2,vde2(iS),veth2,wpi,vith))
        [wrtemp,witemp,kmaxtemp,wrmaxtemp,wimaxtemp,vphmaxtemp,residualtemp] = dg_solver(vde1,vde2(iS),veth1,veth2,vith,wpe,wpe1(iR),wpe2(iR),wpi,kvec,doposwi);
        wr{iR,iS} = wrtemp;
        wi{iR,iS} = witemp;
        residual{iR,iS} = residualtemp;
        wimax(iR,iS) = wimaxtemp;
        wrmax(iR,iS) = wrmaxtemp;
        kmax(iR,iS) = kmaxtemp;
        vphmax(iR,iS) = vphmaxtemp;        
    end
    toc
end

tosave.wr = wr;
tosave.wi = wi;
tosave.residual = residual;
tosave.wimax = wimax;
tosave.wrmax = wrmax;
tosave.kmax = kmax;
tosave.vphmax = vphmax;

plot(hca,kvec*Ld,[witemp;wrtemp]/wpe); legend('w_i','w_r')

drawnow
eval(['tosave' num2str(doposwi) '=tosave'])
end
toc
hold(hca,'off')

%% Save
savePath = '/home/cecilia/Research/BeamSolver/';
saveName = [datestr(now,'yyyy-mm-ddTHHMMSS') '_dg_solver'];
save([savePath saveName])
disp(['saved ' savePath saveName])

%%
if 0
    %%
subplot(2,2,1)
pcolor(R,S*veth2/veth1,wrmax'/wpi)
set(gca,'clim',[0 20]) 
colorbar
subplot(2,2,2)
pcolor(R,S,wimax'/wpi)
set(gca,'clim',[0 5])
colorbar
subplot(2,2,3)
pcolor(R,S,vphmax'/veth1)
set(gca,'clim',[0 0.8])
colorbar
subplot(2,2,4)
pcolor(R,S,kmax'*Ld)
set(gca,'clim',[0 1.2])
colorbar

colormap(cn.cmap('bluered3'))
end
%%
if 0
for ip = 1:nR; h(ip) = subplot(1,nv,ip); end
for ip = 1:nR
   plot(h(ip),kvec*Ld,wr{ip,1},'r',kvec*Ld,wr{ip,2},'b',kvec*Ld,wr{ip,3},'g',...
              kvec*Ld,wi{ip,1},'r--',kvec*Ld,wi{ip,2},'b--',kvec*Ld,wi{ip,3},'g--') 
end
end
%%
if 0
    %%
    iR = 2;
    wr_kS = nan(nk,nv);
    wi_kS = nan(nk,nv);
    for iv = 1:nv;
        wr_kS(1:numel(wr{iR,iv}),iv) = wr{iR,iv};
        wi_kS(1:numel(wr{iR,iv}),iv) = wi{iR,iv};
    end
    
    subplot(2,1,1)
    pcolor(kvec*Ld,S,wr_kS'/wpi)
    set(gca,'clim',[0 30])
    colorbar
    shading flat
    title(['R=' num2str(R(iR))])
    subplot(2,1,2)
    pcolor(kvec*Ld,S,wi_kS'/wpi)
    set(gca,'clim',[0 5])
    colorbar
    shading flat
    
end
%%
if 0
%%

plot(kvec*Ld,wia/wpi,kvec*Ld,wga/wpi)
plot(kvec*Ld,wia/wpe,kvec*Ld,wga/wpe)
xlabel('k\lambda_D')
ylabel('\omega_{max}/\omega_{pi}, \gamma_{max}/\omega_{pi}')
distrStr = {['R=' num2str(ne2/ni,'%.2f')],...
            ['Ti=' num2str(Ti,'%.0f')],...
            ['Te1=' num2str(Te1,'%.0f')],...
            ['Te2=' num2str(Te2,'%.0f')],...
            ['S=' num2str(vd2/veth1,'%.2f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')

%find(min(abs(iadisp(wc))))
end