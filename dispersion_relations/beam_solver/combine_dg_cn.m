if 0 % in case others gets written over
    loadPath = '/Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/';
    matName = 'fig_theory.mat';
end

% load fsolve solver
loadPath = '/Users/Cecilia/Research/EH/BeamSolver/';
matName = '2015-05-17T214209_beam_solver_110kx29Sx20Rx1Te2x1Ti.mat';
load([loadPath matName])

cns.wpe = param.ope*sqrt(1/n); % n=0.06->1
cns.wpi = param.opi*sqrt(1/n); % n=0.06->1
cns.R = R;
cns.S = S;
cns.wi = wimax'*param.wnorm*sqrt(1/n); % (S,R)->(R,S) and n=0.06->1
cns.wr = wrmax'*param.wnorm*sqrt(1/n); % (S,R)->(R,S) and n=0.06->1
cns.vph = vphmax'*param.vte1;
cns.k = kmax';
cns.vtebg = param.vte1;

%% load daniels solver
loadPath = '/Users/Cecilia/Research/EH/BeamSolver/';
matName = '2015-05-13T004113_dg_solver.mat';
load([loadPath matName])

dgs.wpe = wpe*1e3;
dgs.wpi = wpi*1e3;
dgs.R = R(2:end); % take away R~0
dgs.S = S;
dgs.wi = wimax(2:end,:)*dgs.wpi; % take away R~0
dgs.wr = wrmax(2:end,:)*dgs.wpi; % take away R~0
dgs.vph = vphmax(2:end,:); % take away R~0
dgs.k = kmax(2:end,:); % take away R~0
dgs.vtebg = veth1;
%% Combine them
% First, set all damped modes to wi = 0; wr = 0; k = 0; vph = NaN,
s{1} = cns;
s{2} = dgs;
for n = 1:2    
    nanwilim = 1;
    s{n}.wr(s{n}.wi<nanwilim) = NaN;
    s{n}.k(s{n}.wi<nanwilim) = NaN;
    s{n}.vph(s{n}.wi<nanwilim) = NaN;
    s{n}.wi(s{n}.wi<nanwilim) = NaN;  
    s{n}.wi(s{n}.wi>1e6) = NaN;
end

s{3}=s{2};
s0 = nan([size(s{1}.wi) 2]);

% Automatic selection
for iR = 1:numel(s{1}.R);
    for iS = 1:numel(s{1}.S);                    
        [Y,I] = max([s{1}.wi(iR,iS) s{2}.wi(iR,iS)]);
        if iR == 14 && iS == 20; I = 1; end %  avvikande punkten iR=14 ist. f?r 15 d? jag tog bort R=1 indexet fr?n dg.       
        if s{1}.wi(iR,iS)<100 && s{1}.wr(iR,iS)>1e5; I = 2; end
        %disp([num2str(I) ' ' num2str(Y)])
        s{3}.wi(iR,iS) = s{I}.wi(iR,iS);
        s{3}.wr(iR,iS) = s{I}.wr(iR,iS);
        s{3}.k(iR,iS) = s{I}.k(iR,iS);
        s{3}.vph(iR,iS) = s{I}.vph(iR,iS);
        s{3}.ind(iR,iS) = I;               
    end
end
for ii=1:3;
    s{ii}.wi(1,12)=NaN; s{ii}.wr(1,12)=NaN; s{ii}.k(1,12)=NaN; s{ii}.vph(1,12)=NaN;
end


% Manual selection
s{4}=s{2};
s{4}.ind = 2*ones(numel(s{1}.R),numel(s{1}.S));  
if 0
s{4}.ind = 2*ones(numel(s{1}.R),numel(s{1}.S));  
s{4}.R = s{1}.R; s{4}.S = s{1}.S;
s{4}.wpi = s{1}.wpi; s{4}.wpe = s{1}.wpe;
s{4}.vtebg = s{1}.vtebg;
s{4}.vph = s{1}.vph
end
Ss = [12 8  10 12 14 15 14 11 12 23 13 12 10 10 11 12 16 19 27 20 17 15 14 15 16 17]; % change these indices to ones
Rs = [13 14 14 15 15 15 16 16 17 18 18 18 18 19 19 19 19 19 19 14 11 18 17 17 17 17]; % change these indices to ones

for kk = 1:numel(Ss); 
    s{4}.ind(Rs(kk),Ss(kk)) = 1; 
    s{4}.wi(Rs(kk),Ss(kk)) = s{1}.wi(Rs(kk),Ss(kk));
    s{4}.wr(Rs(kk),Ss(kk)) = s{1}.wr(Rs(kk),Ss(kk));
    s{4}.k(Rs(kk),Ss(kk)) = s{1}.k(Rs(kk),Ss(kk))/(Ld);
    s{4}.vph(Rs(kk),Ss(kk)) = s{1}.vph(Rs(kk),Ss(kk));    
end

if 0
for iR = 1:numel(s{1}.R);
    for iS = 1:numel(s{1}.S);                    
        %disp([num2str(I) ' ' num2str(Y)])
        s{4}.wi(iR,iS) = s{s{4}.ind(iR,iS)}.wi(iR,iS);
        s{4}.wr(iR,iS) = s{s{4}.ind(iR,iS)}.wr(iR,iS);
        %s{4}.k(iR,iS) = s{s{4}.ind(iR,iS)}.k(iR,iS);     
        s{s{4}.ind(iR,iS)}.k(iR,iS)
        %s{4}.vph(iR,iS) = s{s{4}.ind(iR,iS)}.wr(iR,iS)./s{s{4}.ind(iR,iS)}.k(iR,iS);
        s{4}.vph(iR,iS) = s{s{4}.ind(iR,iS)}.vph(iR,iS);
    end
end
end
%s{4}.vph = s{4}.wr./s{4}.k;
% Make some values NaN
Rnan = [1 1 1 1 1 1  1  1  1  1  1  1  1  2];
Snan = [3 4 7 8 9 10 11 12 13 14 15 16 17 9];
for kk = 1:numel(Snan); 
    s{4}.vph(Rnan(kk),Snan(kk)) = NaN; 
    s{4}.wr(Rnan(kk),Snan(kk)) = NaN; 
    s{4}.wi(Rnan(kk),Snan(kk)) = NaN; 
    s{4}.k(Rnan(kk),Snan(kk)) = NaN; 
end

s{4}.wi(1,12)=0; s{4}.wr(1,12)=0; s{4}.k(1,12)=0; s{4}.vph(1,12)=NaN;
for ii = 3:4; for iS=2:3 for iR = 1:numel(s{1}.R);    
    s{ii}.wi(iR,iS)=NaN; s{ii}.wr(iR,iS)=NaN; s{ii}.k(iR,iS)=NaN; s{ii}.vph(iR,iS)=NaN;
end; end; end
iR=5; iS = 5;
s{4}.wi(iR,iS)=NaN; s{4}.wr(iR,iS)=NaN; s{4}.k(iR,iS)=NaN; s{4}.vph(iR,iS)=NaN;
iR=1; iS = 12;
s{4}.wi(iR,iS)=NaN; s{4}.wr(iR,iS)=NaN; s{4}.k(iR,iS)=NaN; s{4}.vph(iR,iS)=NaN;
iR=1; iS = 5;
s{4}.wi(iR,iS)=NaN; s{4}.wr(iR,iS)=NaN; s{4}.k(iR,iS)=NaN; s{4}.vph(iR,iS)=NaN;

% kS of wr wi

iR = 7; iRplot=iR;
wr_kS = nan(nk,nv);
wi_kS = nan(nk,nv);
for iS = 1:nv;
    wr_kS(1:numel(wr{iR,iS}),iS) = wr{iR,iS};
    wi_kS(1:numel(wr{iR,iS}),iS) = wi{iR,iS};
end
if iRplot == 5
iS = 2;
nkk = numel(wr{iR,iS});
wr_kS(1:nkk,iS) = 0.5*(wr{iR,iS+1}(1:nkk)+wr{iR,iS-1}(1:nkk))';
wi_kS(1:nkk,iS) = 0.5*(wi{iR,iS+1}(1:nkk)+wi{iR,iS-1}(1:nkk))';
end
if iRplot == 6
iS = 2;
nkk = numel(wr{iR,iS});
wr_kS(1:nkk,iS) = 0.5*(wr{iR,iS+1}(1:nkk)+wr{iR,iS-1}(1:nkk))';
wi_kS(1:nkk,iS) = 0.5*(wi{iR,iS+1}(1:nkk)+wi{iR,iS-1}(1:nkk))';
iS = 5;
nkk = numel(wr{iR,iS}); nkk = 200;
wr_kS(1:nkk,iS) = 0.5*(wr{iR,iS+1}(1:nkk)+wr{iR,iS-1}(1:nkk))';
wi_kS(1:nkk,iS) = 0.5*(wi{iR,iS+1}(1:nkk)+wi{iR,iS-1}(1:nkk))';
iS = 2;
nkk = numel(wr{iR,iS}); nkk = 200;
wr_kS(1:nkk,iS) = 0.54*(wr{iR,4}(1:nkk)+wr{iR,1}(1:nkk))';
wi_kS(1:nkk,iS) = 0.54*(wi{iR,4}(1:nkk)+wi{iR,1}(1:nkk))';
iS = 3;
nkk = numel(wr{iR,iS}); nkk = 200;
wr_kS(1:nkk,iS) = 0.47*(wr{iR,4}(1:nkk)+wr{iR,1}(1:nkk))';
wi_kS(1:nkk,iS) = 0.47*(wi{iR,4}(1:nkk)+wi{iR,1}(1:nkk))';
end
if iRplot == 7
    iS = 2;
    nkk = numel(wr{iR,iS});
    wr_kS(1:nkk,iS) = 0.5*(wr{iR,iS+1}(1:nkk)+wr{iR,iS-1}(1:nkk))';
    wi_kS(1:nkk,iS) = 0.5*(wi{iR,iS+1}(1:nkk)+wi{iR,iS-1}(1:nkk))';
end
%% Make figure
ii = 4; % 1 = cn, 2 = dg, 3 = automatic, 4 = manual
set(0,'defaultAxesFontSize',18);
set(0,'DefaultTextFontSize',18);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');
nPlots = 6;
nrows = 3;
for kk = 1:nPlots; h(kk)=subplot(nrows,2,kk); end
colormap(cn.cmap('bluered3'))
isub=1;    
% figure limits
limk = [0 1];
limwikS = [-2 5];
limwrkS = [0 60]; limwrKSlog = [];
limwi = [0 7]; 
limwr = [0 30];
limv = [0 1.2];
limR = [0.05 0.99];
limS = [0.1 1.5];
liml = [0 20];

doLog = 1;
Rg = [0 0.075:0.05:1 1];
Sg = 0.075:0.05:1.55;
if 1 % wi vs k,S 
    hca = h(isub); isub=isub+1;           
    pcolor(hca,kvec*Ld,S,wi_kS'/wpi)
    xlabel(hca,'k\lambda_{De}'); ylabel(hca,'S')    
    ch = colorbar('peer',hca);                
    ylabel(ch,'\omega_{i}/\omega_{pi}')      
    set(hca,'clim',limwikS(2)*[-1 1],'ylim',limS,'xlim',limk); 
    set(ch,'ylim',limwikS)
    shading(hca,'flat')    
    xlab{isub-1} = [' R = ' num2str(R(iRplot))];
end
if 1 % wr vs k,S 
    if doLog
        hca = h(isub); isub=isub+1;           
        pcolor(hca,kvec*Ld,S,log10(wr_kS'/wpi))
        xlabel(hca,'k\lambda_{De}'); ylabel(hca,'S') 
        ch = colorbar('peer',hca);             
        ylabel(ch,'\omega_{r}/\omega_{pi}')
        set(hca,'clim',2*[-1 1],'ylim',limS,'xlim',limk)
        set(ch,'ylim',2*[-1 1])
        ctick = -2:2; cticklabelnum = 10.^ctick;        
        %for kk = 1:numel(ctick)
        %    cticklabels{kk} = ['10^' num2str(ctick(kk))];
        %end
        cticklabels = [0.01 0.1 1 10 100];
        set(ch,'ytick',ctick,'yticklabel',cticklabels)
        shading(hca,'flat')   
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
    else
        hca = h(isub); isub=isub+1;           
        pcolor(hca,kvec*Ld,S,wr_kS'/wpi)
        xlabel(hca,'k\lambda_{De}'); ylabel(hca,'S') 
        ch = colorbar('peer',hca);             
        ylabel(ch,'\omega_{r}/\omega_{pi}')
        set(hca,'clim',limwrkS(2)*[-1 1],'ylim',limS,'xlim',limk)
        set(ch,'ylim',limwrkS)
        shading(hca,'flat')   
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
    end
end

if 1 % wimax vs R,S, log, growthrate as colorscale
    hca = h(isub); isub=isub+1;         
    pcolor(hca,s{ii}.R,s{ii}.S,s{ii}.wi'/s{ii}.wpi)
    %surf(hca,Rg,Sg,s{ii}.wi'/s{ii}.wpi)               
    xlabel(hca,'R'); ylabel(hca,'S'); %title(hca,'Maximum growth rate')
    ch = colorbar('peer',hca);  ylabel(ch,'\omega_{i}/\omega_{pi}')  
    shading(hca,'flat'); box(hca,'on')
    set(hca,'clim',limwi(2)*[-1 1],'ylim',limS,'xlim',limR)
    set(ch,'ylim',limwi)    
end
if 1 % wrmax vs R,S, log, growthrate as colorscale
    if doLog
        hca = h(isub); isub=isub+1;             
        pcolor(hca,s{ii}.R,s{ii}.S,log10(s{ii}.wr'/s{ii}.wpi))
        %pcolor(hca,Rg,Sg,s{ii}.wr'/s{ii}.wpi)               
        xlabel(hca,'R'); ylabel(hca,'S')
        ch = colorbar('peer',hca); ylabel(ch,'\omega_{r}/\omega_{pi}')
        %set(ch,'cscale','log')
        shading(hca,'flat'); box(hca,'on')
        cticklabels = [ 0.1 1 10 100];
        ctick = -2:2; cticklabelnum = 10.^ctick;   
        %for kk = 1:numel(ctick)
        %    texcticklabels{kk} = ['10^' num2str(ctick(kk))];
        %end
        %set(ch,'TickLabelInterpreter','tex','yticklabel',texcticklabels)
        set(ch,'ytick',ctick,'yticklabel',cticklabels)
        %ctick
        chh=ch;
        set(hca,'clim',[-1 2],'ylim',limS,'xlim',limR)
        set(ch,'ylim',[-1 2])        
    else
        hca = h(isub); isub=isub+1;             
        pcolor(hca,s{ii}.R,s{ii}.S,s{ii}.wr'/s{ii}.wpi)               
        %pcolor(hca,Rg,Sg,s{ii}.wr'/s{ii}.wpi)               
        xlabel(hca,'R'); ylabel(hca,'S')
        ch = colorbar('peer',hca); ylabel(ch,'\omega_{r}/\omega_{pi}')
        %set(ch,'cscale','log')
        shading(hca,'flat'); box(hca,'on')
        set(hca,'clim',limwr(2)*[-1 1],'ylim',limS,'xlim',limR)
        set(ch,'ylim',limwr)        
    end
end
if 1 % vph vs R,S, log, growthrate as colorscale
    if 0
        hca = h(isub); isub=isub+1;     
        pcolor(hca,s{ii}.R,s{ii}.S,log10(s{ii}.vph'/s{4}.vtebg))
        %pcolor(hca,Rg,Sg,s{ii}.vph'/s{4}.vtebg) 
        xlabel(hca,'R'); ylabel(hca,'S')
        set(hca,'ylim',limS,'xlim',limR,'clim',[-2 0.2],'zscale','lin') 
        caxis(hca,10.^[-2 0.2]);
        ch = colorbar('peer',hca,'yscale','log');    
        cticklabels = [0.001 0.01 0.1 1 10 100];
        ctick = log10(cticklabels); cticklabelnum = 10.^ctick;   
        set(ch,'ylim',10.^[-2 0.2],'ytick',[10.^[-2 -1 1]])
        ylabel(ch,'v_{ph}/v_{te,bg}')                        
        shading(hca,'flat'); box(hca,'on');
        caxis(hca,10.^[-2 0.2]);
    elseif doLog
        hca = h(isub); isub=isub+1;     
        pcolor(hca,s{ii}.R,s{ii}.S,log10(s{ii}.vph'/s{4}.vtebg))
        %pcolor(hca,Rg,Sg,s{ii}.vph'/s{4}.vtebg) 
        xlabel(hca,'R'); ylabel(hca,'S')
        set(hca,'ylim',limS,'xlim',limR,'clim',[-2.2 0.2],'zscale','lin')
        %set(hca,'zlim',[-0.01 0.12])
        ch = colorbar('peer',hca);    
        cticklabels = [0.001 0.01 0.1 1 10 100];
        ctick = log10(cticklabels); cticklabelnum = 10.^ctick;   
        set(ch,'ylim',[-2.2 0.2],'ytick',ctick,'yticklabel',cticklabels)
        ylabel(ch,'v_{ph}/v_{te,bg}')                        
        shading(hca,'flat'); box(hca,'on');
    else
        hca = h(isub); isub=isub+1;     
        pcolor(hca,s{ii}.R,s{ii}.S,s{ii}.vph'/s{4}.vtebg) 
        %pcolor(hca,Rg,Sg,s{ii}.vph'/s{4}.vtebg) 
        xlabel(hca,'R'); ylabel(hca,'S')
        set(hca,'ylim',limS,'xlim',limR,'clim',limv(2)*[-1 1],'zscale','lin')
        %set(hca,'zlim',[-0.01 0.12])
        ch = colorbar('peer',hca);    
        set(ch,'ylim',limv)                        
        ylabel(ch,'v_{ph}/v_{te,bg}')                
        %set(ch,'ylim',clim)
        shading(hca,'flat'); box(hca,'on');
    end
end
if 1 % lambda vs R,S, log, growthrate as colorscale
    hca = h(isub); isub=isub+1;     
    %pcolor(hca,s{ii}.R,s{ii}.S,2*pi./s{ii}.k') 
    kLd =  s{4}.k*Ld; % gives the k*lambda_De seen in panel (a) for ex.
    LLd = 2*pi./kLd;
    lambda_De = s{ii}.vtebg/sqrt(2)/s{ii}.wpe;
    Lplot = 2*pi./s{ii}.k'/1e3;
    pcolor(hca,s{ii}.R,s{ii}.S,LLd') 
    %pcolor(hca,Rg,Sg,2*pi./s{ii}.k'/s{ii}.vtebg/sqrt(2)*s{ii}.wpe*1e-3) 
    xlabel(hca,'R'); ylabel(hca,'S')
    set(hca,'ylim',limS,'xlim',limR,'clim',liml(2)*[-1 1],'zscale','lin')
    %set(hca,'zlim',[-0.01 0.12])
    ch = colorbar('peer',hca);    
    set(ch,'ylim',liml)
    ylabel(ch,'\lambda\omega_{pe}/v_{te,bg}{2}^{1/2}'); ylabel(ch,'\lambda/\lambda_{De}')                
    %set(ch,'ylim',clim)
    shading(hca,'flat'); box(hca,'on');
end
if 0
    hca = h(isub); isub=isub+1;             
    pcolor(hca,s{ii}.R,s{ii}.S,s{ii}.ind'); 
    ch = colorbar('peer',hca);
end
if 0 % wi vs k,S 
    hca = h(isub); isub=isub+1;   
    RSratio=R'*(1./(S));
    RSratio=R'*(1./(S-0.3));
    pcolor(hca,R,S,(RSratio'))
    xlabel(hca,'R'); ylabel(hca,'S')    
    ch = colorbar('peer',hca);                
    ylabel(ch,'R/S')      
    set(hca,'clim',1+1*[-1 1],'ylim',limS,'xlim',limR); 
    set(ch,'ylim',1+1*[-1 1])
    shading(hca,'flat')    
    %xlab{isub-1} = [' R = ' num2str(R(iRplot))];
end

if 0 % wi vs k,S 
    hca = h(isub); isub=isub+1;   
    SRratio=(1./R)'*S/1.5;
    SRratio=(1./R)'*(S-0.3);
    pcolor(hca,R,S,(SRratio'))
    xlabel(hca,'R'); ylabel(hca,'S')    
    ch = colorbar('peer',hca);                
    ylabel(ch,'S/R')      
    set(hca,'clim',1+3*[-1 1],'ylim',limS,'xlim',limk); 
    set(ch,'ylim',1+3*[-1 1])
    shading(hca,'flat')    
    %xlab{isub-1} = [' R = ' num2str(R(iRplot))];
end
%colormap(cn.cmap('whitered'));

    % add labels
abc = {'a)','b)','c)','d)','e)','f)','g)','h)'};
for kk = 1:nPlots        
    xl = get(h(kk),'xlim');
    yl = get(h(kk),'ylim');
    axes(h(kk))
    try abc{kk} = [abc{kk} xlab{kk}]; end
    %text(xl(1)+diff(xl)*0.05,yl(2)-diff(xl)*0.05,abc{kk},'horizontalalignment','left','verticalalignment','top','color','white')
    %text(xl(2)-diff(xl)*0.05,yl(2)-diff(xl)*0.05,abc{kk},'horizontalalignment','right','verticalalignment','top','color','white')
    text(xl(1)+diff(xl)*0.03,yl(1)+diff(xl)*0.03,abc{kk},'horizontalalignment','left','verticalalignment','bottom')
    %text(xl(2)-diff(xl)*0.03,yl(1)+diff(xl)*0.03,abc{kk},'horizontalalignment','right','verticalalignment','bottom')
    axis(h(kk),'square')
end


% change position of axes
for kk = 1:2:nPlots    
    pos = get(h(kk),'position');
    set(h(kk),'position',[pos(1)-0.03 pos(2) pos(3)*1.1 pos(4)*1.1])    
end
for kk = 2:2:nPlots
    pos = get(h(kk),'position');
    set(h(kk),'position',[pos(1) pos(2) pos(3)*1.1 pos(4)*1.1])    
end

%%
subplot(2,2,1)
pcolor(s{4}.R,s{4}.S,s{4}.vph'./vT')
cb=colorbar;
set(cb,'ytick',-2:0.2:2)
ylabel(cb,'v_{ph}/v_{T}')
set(gca,'clim',1+1*[-1 1])
xlabel(gca,'R')
ylabel(gca,'S')
axis square

subplot(2,2,2)
pcolor(s{4}.R,s{4}.S,vT'./s{4}.vph')
cb=colorbar;
ylabel(cb,'v_T/v_{ph}')
set(gca,'clim',0+27*[0 1])
xlabel(gca,'R')
ylabel(gca,'S')
axis square

subplot(2,2,3)
pcolor(s{4}.R,s{4}.S,vT'./veth2)
cb=colorbar;
ylabel(cb,'v_T/v_{te,beam}')
set(gca,'clim',0+7*[0 1])
xlabel(gca,'R')
ylabel(gca,'S')
axis square

colormap(cn.cmap('bluered3'))
%%
wrlim = [0 7];
subplot(4,1,1)
ind=1;
pcolor(s{ind}.R,s{ind}.S,s{ind}.wi'/s{ind}.wpi); hold on;
%surf(dgs.R,dgs.S,dgs.wi'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off

subplot(4,1,2)
ind=2;
pcolor(s{ind}.R,s{ind}.S,s{ind}.wi'/s{ind}.wpi); hold on;
%surf(dgs.R,dgs.S,dgs.wi'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off

subplot(4,1,3)
ind=3;
pcolor(s{ind}.R,s{ind}.S,s{ind}.wi'/s{ind}.wpi); hold on;
%surf(dgs.R,dgs.S,dgs.wi'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off

subplot(4,1,4)
ind=3;
pcolor(s{1}.R,s{1}.S,s{1}.wi'/s{1}.wpi-s{2}.wi'/s{2}.wpi); hold on;
%surf(dgs.R,dgs.S,dgs.wi'); hold on;
colorbar
set(gca,'clim',0.01*[-1 1])
hold off
%%
subplot(3,1,2)
%surf(cns.R,cns.S,cns.wi); hold on;
pcolor(dgs.R(1:end),dgs.S,dgs.wr(1:end,:)'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off
subplot(3,1,3)
%surf(cns.R,cns.S,cns.wi); hold on;
surf(dgs.R(2:end),dgs.S,dgs.wr(2:end,:)'); hold on;
surf(cns.R,cns.S,cns.wr/cns.wpi); hold on;
colorbar
set(gca,'clim',wrlim)
hold off
%%
wrlim = [0 35];
subplot(3,1,1)
pcolor(cns.R,cns.S,cns.wr/cns.wpi); hold on;
%surf(dgs.R,dgs.S,dgs.wi'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off
subplot(3,1,2)
%surf(cns.R,cns.S,cns.wi); hold on;
pcolor(dgs.R(1:end),dgs.S,dgs.wr(1:end,:)'); hold on;
colorbar
set(gca,'clim',wrlim)
hold off
subplot(3,1,3)
%surf(cns.R,cns.S,cns.wi); hold on;
surf(dgs.R(2:end),dgs.S,dgs.wr(2:end,:)'); hold on;
surf(cns.R,cns.S,cns.wr/cns.wpi); hold on;
colorbar
set(gca,'clim',wrlim)
hold off
%%
s{4}.vT = repmat(s{4}.S,numel(s{4}.R),1)*s{4}.vtebg-s{4}.vph;

px = s{4}.vph./s{4}.vtebg);
py = s{4}.vph./s{4}.vT;
plot(,'o')

%%
nPlots=2;
nrows = 1;
for kk = 1:nPlots; h(kk)=subplot(nrows,2,kk); end
isub = 1;
if 1 % wi vs k,S 
    hca = h(isub); isub=isub+1;   
    RSratio=R'*(1./S);
    pcolor(hca,R,S,log10(RSratio'))
    xlabel(hca,'R'); ylabel(hca,'S')    
    ch = colorbar('peer',hca);                
    ylabel(ch,'R/S')      
    %set(hca,'clim',limwikS(2)*[-1 1],'ylim',limS,'xlim',limk); 
    %set(ch,'ylim',limwikS)
    shading(hca,'flat')    
    %xlab{isub-1} = [' R = ' num2str(R(iRplot))];
end

if 1 % wi vs k,S 
    hca = h(isub); isub=isub+1;   
    SRratio=(1./R)'*S;
    pcolor(hca,R,S,log10(SRratio'))
    xlabel(hca,'R'); ylabel(hca,'S')    
    ch = colorbar('peer',hca);                
    ylabel(ch,'S/R')      
    %set(hca,'clim',limwikS(2)*[-1 1],'ylim',limS,'xlim',limk); 
    %set(ch,'ylim',limwikS)
    shading(hca,'flat')    
    %xlab{isub-1} = [' R = ' num2str(R(iRplot))];
end
