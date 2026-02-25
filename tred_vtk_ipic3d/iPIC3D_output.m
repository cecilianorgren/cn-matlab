%filepath = '/Users/Cecilia/iPIC3D/hdf5/';
filepath = '/Users/Cecilia/iPIC3D/newiPIC3D/ipic3d-kth/iPic3D-KTH/build/';
%subfolder = 'old_data/ampere_test_2/';
subfolder = 'data/';
subfolder = 'old_data/debugged_code_ivy/';
filename = 'proc0.hdf';
%h5disp([filepath filename],'/particles/species_2/','min');

% get cycles and dt
info = h5info([filepath subfolder 'proc0.hdf'],'/moments/species_2/Jx');
ncycles = size(info.Datasets,1);
%[a{1:ncycles,1}] = deal(info.Datasets.Name);
cycles = zeros(ncycles,1);
for kk = 1:ncycles
    %info.Datasets(kk).Name(7:end);
    cycles(kk) = str2double(info.Datasets(kk).Name(7:end));
end
cycles = sort(cycles);

fid = fopen([filepath subfolder 'SimulationData.txt']);
simdata = fscanf(fid,'%s',inf);
ind1 = strfind(simdata,'Timestep=')+9;
ind2 = strfind(simdata,'Numberofcycles=')-1;
dt = str2double(simdata(ind1:ind2));
t = dt*cycles;

nspecies = str2double(simdata(strfind(simdata,'Numberofspecies=')+16));

for kk = 1:nspecies
    eval(['ind_a = strfind(simdata,''rhoinitspecies' num2str(kk-1) '='');']);
    eval(['ind_b = strfind(simdata,''rhoinjectspecies' num2str(kk-1) '='');']);
    eval(['n'  num2str(kk-1) '= str2double(simdata(ind_a+16:ind_b-1))'])   
    %eval(['ind_c = strfind(simdata,''qom[%d]=' num2str(kk-1) '='');']);
    %eval(['ind_d = strfind(simdata,''rhoinitspecies' num2str(kk-1) '='');']);
end
%% normalization
[expr,inds] = regexp(simdata,'qom','match');
inds_a = regexp(simdata(inds(2)+[8:17]),'\d');

units = irf_units;
%%

% load data

% More easily adaptable
if 0
%%
%cycles = [0; 100; 200; 300; 400; 500; 600; 700; 800; 900];
%cycles = 0:2:28;
ncycles = numel(cycles);
nSpecies = 4;
nPanels = nSpecies+1;
colors = [0 1 0;1 0 0;0 0 1];
for kk=1:nPanels; h(kk) = subplot(1,nPanels,kk); end
for ic = 1:ncycles
    for is = 0:(nSpecies-1)
        eval(['u' num2str(is) '=h5read([ filepath subfolder filename ],''/particles/species_' num2str(is) '/u/cycle_' num2str(cycles(ic)) ''');'])
        eval(['x' num2str(is) '=h5read([ filepath subfolder filename ],''/particles/species_' num2str(is) '/x/cycle_' num2str(cycles(ic)) ''');'])
        eval(['u' num2str(is) '=h5read([ filepath subfolder filename ],''/particles/species_' num2str(is) '/u/cycle_' num2str(cycles(ic)) ''');'])
    end
    
    isub=1;
    for is = 0:(nSpecies-1)
        hca=h(isub); isub=isub+1;
        eval(['x_tp = x' num2str(is) ';']);
        eval(['u_tp = u' num2str(is) ';']);
        plot(hca,x_tp,u_tp,'.','color',colors(is+1,:))
        set(hca,'ylim',0.3*[-1 1],'xlim',[0 0.257]);
    end
    
    if 0
    title(hca,'species\_0')   
    plot(h(1),x0,u0,'g.',[0 0.4],mean(u0)*[1 1],'k-','linewidth',2);    
    title(h(2),'species\_1')
    plot(h(2),x1,u1,'r.',[0 0.4],mean(u1)*[1 1],'k-','linewidth',2);
    title(h(3),'species\_2')
    plot(h(3),x2,u2,'b.',[0 0.4],mean(u2)*[1 1],'k-','linewidth',2); 
    plot(h(4),[x0; x2],[u0; u2],'g.')
    plot(h(4),x2,u2,'b.',x1,u1,'r.',x0,u0,'g.')
    plot(h(5),x2,u2-mean(u0),'b.',x1,u1,'r.',x0,u0-mean(u0),'g.')
    for kk=1:5; set(h(kk),'ylim',0.2*[-1 1],'xlim',[0 0.257]); end
    end
    drawnow
    pause(0.8)
end
end
if 0
%%
%cycles = [0; 100; 200; 300; 400; 500; 600; 700; 800; 900];
%cycles = 0:2:28;
%cycles= 0:50:400;
nc = numel(cycles);

% Plot x vx phasespace
nSpecies = 3;
nPanels=5;
for kk=1:nPanels; h(kk) = subplot(1,nPanels,kk); end
for ic = 1:nc
    u0 = h5read([filepath subfolder filename],['/particles/species_0/u/cycle_' num2str(cycles(ic))]);
    v0 = h5read([filepath subfolder filename],['/particles/species_0/v/cycle_' num2str(cycles(ic))]);
    w0 = h5read([filepath subfolder filename],['/particles/species_0/w/cycle_' num2str(cycles(ic))]);
    x0 = h5read([filepath subfolder filename],['/particles/species_0/x/cycle_' num2str(cycles(ic))]);
    y0 = h5read([filepath subfolder filename],['/particles/species_0/y/cycle_' num2str(cycles(ic))]);
    z0 = h5read([filepath subfolder filename],['/particles/species_0/z/cycle_' num2str(cycles(ic))]);
    u1 = h5read([filepath subfolder filename],['/particles/species_1/u/cycle_' num2str(cycles(ic))]);
    v1 = h5read([filepath subfolder filename],['/particles/species_1/v/cycle_' num2str(cycles(ic))]);
    w1 = h5read([filepath subfolder filename],['/particles/species_1/w/cycle_' num2str(cycles(ic))]);
    x1 = h5read([filepath subfolder filename],['/particles/species_1/x/cycle_' num2str(cycles(ic))]);
    y1 = h5read([filepath subfolder filename],['/particles/species_1/y/cycle_' num2str(cycles(ic))]);
    z1 = h5read([filepath subfolder filename],['/particles/species_1/z/cycle_' num2str(cycles(ic))]);
    u2 = h5read([filepath subfolder filename],['/particles/species_2/u/cycle_' num2str(cycles(ic))]);
    v2 = h5read([filepath subfolder filename],['/particles/species_2/v/cycle_' num2str(cycles(ic))]);
    w2 = h5read([filepath subfolder filename],['/particles/species_2/w/cycle_' num2str(cycles(ic))]);
    x2 = h5read([filepath subfolder filename],['/particles/species_2/x/cycle_' num2str(cycles(ic))]);
    y2 = h5read([filepath subfolder filename],['/particles/species_2/y/cycle_' num2str(cycles(ic))]);
    z2 = h5read([filepath subfolder filename],['/particles/species_2/z/cycle_' num2str(cycles(ic))]);
    
    title(h(1),'species\_0')   
    plot(h(1),x0,u0,'g.',[0 0.4],mean(u0)*[1 1],'k-','linewidth',2);    
    title(h(2),'species\_1')
    plot(h(2),x1,u1,'r.',[0 0.4],mean(u1)*[1 1],'k-','linewidth',2);
    title(h(3),'species\_2')
    plot(h(3),x2,u2,'b.',[0 0.4],mean(u2)*[1 1],'k-','linewidth',2); 
    plot(h(4),[x0; x2],[u0; u2],'g.')
    plot(h(4),x2,u2,'b.',x1,u1,'r.',x0,u0,'g.')
    plot(h(5),x2,u2-mean(u0),'b.',x1,u1,'r.',x0,u0-mean(u0),'g.')
    for kk=1:5; set(h(kk),'ylim',0.2*[-1 1],'xlim',[0 0.257]); end
    drawnow
    pause(0.8)
end

end
%% Plot Ex
nPlots=3;
for kk=1:nPlots; h(kk) = subplot(nPlots,1,kk); end
for ic = 1:ncycles   
    tnow = dt*cycles(ic)
    Ex = h5read([filepath subfolder filename],['/fields/Ex/cycle_' num2str(cycles(ic))]);
    Ey = h5read([filepath subfolder filename],['/fields/Ey/cycle_' num2str(cycles(ic))]);
    Ez = h5read([filepath subfolder filename],['/fields/Ez/cycle_' num2str(cycles(ic))]);
    u0 = h5read([filepath subfolder filename],['/particles/species_0/u/cycle_' num2str(cycles(ic))]);
    x0 = h5read([filepath subfolder filename],['/particles/species_0/x/cycle_' num2str(cycles(ic))]);
    u1 = h5read([filepath subfolder filename],['/particles/species_1/u/cycle_' num2str(cycles(ic))]);
    x1 = h5read([filepath subfolder filename],['/particles/species_1/x/cycle_' num2str(cycles(ic))]);
    u2 = h5read([filepath subfolder filename],['/particles/species_2/u/cycle_' num2str(cycles(ic))]);
    x2 = h5read([filepath subfolder filename],['/particles/species_2/x/cycle_' num2str(cycles(ic))]);
    % density
    n0 = h5read([filepath subfolder filename],['/moments/species_0/rho/cycle_' num2str(cycles(ic))]);
    n1 = h5read([filepath subfolder filename],['/moments/species_1/rho/cycle_' num2str(cycles(ic))]);
    n2 = h5read([filepath subfolder filename],['/moments/species_2/rho/cycle_' num2str(cycles(ic))]);
    
    xlim = [0 0.257];
    xs = linspace(xlim(1),xlim(2),size(n0,3));
    isub=1;
    if 0
        hca = h(isub); isub=isub+1;
        pcolor(hca,squeeze(Ex(1,:,:))); 
        title(hca,'Ex'); 
        ch1=colorbar('peer',hca);
        shading(hca,'flat')
        set(hca,'clim',0.02*[-1 1]);
        xlabel(hca,'x')
    end
    if 1
        hca = h(isub); isub=isub+1;
        plot(hca,xs,squeeze(Ex(1,1,:)));
        %title(hca,['t = ' num2str(tnow) ' \omega_{pi} = ' num2str(t*sqrt(units.mp/units.me*R),'%.1f') ' \omega_{pe}']); 
        set(hca,'ylim',0.3*[-1 1],'xlim',xlim);
        xlabel(hca,'x')
        xlabel(hca,'E_x')
    end
    if 1
        hca = h(isub); isub=isub+1;
        plot(hca,xs,squeeze(n0(1,1,:)),...
                 xs,squeeze(n1(1,1,:)),...
                 xs,squeeze(n2(1,1,:)),...
                 xs,squeeze(n0(1,1,:)+n2(1,1,:)),...
                 xlim,mean(squeeze(n0(1,1,:)+n1(1,1,:)+n2(1,1,:)))*[1 1],'k');
       % legend(hca,'e_{bg}','i','e_{beam}','e_{tot}','e+i','location','bestoutside')
        ylabel(hca,'n'); 
        set(hca,'ylim',0.1*[-1 1],'xlim',xlim);
        xlabel(hca,'x')
    end
    if 1
        hca = h(isub); isub=isub+1;
        plot(hca,x0,u0,'g.',x1,u1,'r.',x2,u2,'b.')    
        set(hca,'ylim',0.2*[-1 1],'xlim',xlim);    
        xlabel(hca,'x')
        ylabel('f(v_x)')
    end
    drawnow
    pause(0.5)
end
%% compare the phases of E and ve
mEx = zeros(ncycles,1);
mEy = zeros(ncycles,1);
mEz = zeros(ncycles,1);
mBx = zeros(ncycles,1);
mBy = zeros(ncycles,1);
mBz = zeros(ncycles,1);
mJx0 = zeros(ncycles,1);
mJx1 = zeros(ncycles,1);
mJx2 = zeros(ncycles,1);
mu0 = zeros(ncycles,1);
mu1 = zeros(ncycles,1);
mu2 = zeros(ncycles,1);


for ic = 1:ncycles
    Ex = h5read([filepath subfolder filename],['/fields/Ex/cycle_' num2str(cycles(ic))]);        
    mEx(ic) = mean(Ex(1,1,:));
    Bx = h5read([filepath subfolder filename],['/fields/Bx/cycle_' num2str(cycles(ic))]);        
    mBy(ic) = mean(Bx(1,1,:));
    By = h5read([filepath subfolder filename],['/fields/By/cycle_' num2str(cycles(ic))]);        
    mBy(ic) = mean(By(1,1,:));
    Bz = h5read([filepath subfolder filename],['/fields/Bz/cycle_' num2str(cycles(ic))]);        
    mBz(ic) = mean(Bz(1,1,:));
    try
    Jx3 = h5read([filepath subfolder filename],['/moments/species_3/Jx/cycle_' num2str(cycles(ic))]);
    mJx3(ic) = mean(Jx3(1,1,:));
    mu3(ic) = mean(h5read([filepath subfolder filename],['/particles/species_3/u/cycle_' num2str(cycles(ic))]));    
    end
    Jx2 = h5read([filepath subfolder filename],['/moments/species_2/Jx/cycle_' num2str(cycles(ic))]);
    mJx2(ic) = mean(Jx2(1,1,:));
    Jx1 = h5read([filepath subfolder filename],['/moments/species_1/Jx/cycle_' num2str(cycles(ic))]);
    mJx1(ic) = mean(Jx1(1,1,:));
    Jx0 = h5read([filepath subfolder filename],['/moments/species_0/Jx/cycle_' num2str(cycles(ic))]);
    mJx0(ic) = mean(Jx0(1,1,:));
    mu0(ic) = mean(h5read([filepath subfolder filename],['/particles/species_0/u/cycle_' num2str(cycles(ic))]));    
    mu1(ic) = mean(h5read([filepath subfolder filename],['/particles/species_1/u/cycle_' num2str(cycles(ic))]));    
    mu2(ic) = mean(h5read([filepath subfolder filename],['/particles/species_2/u/cycle_' num2str(cycles(ic))]));    
end

%% Get time derivative of Ex

x = tocolumn(t);
y = tocolumn(mEx);
type = fittype('a*sin(x*omega)');
Efit=fit(tocolumn(t),tocolumn(mEx),type,'Start',[1 1]);
coeffs = coeffvalues(Efit); E0 = coeffs(1); omega = coeffs(2);
plot(x, y, 'r*', linspace(t(1),t(end),100), Efit(linspace(t(1),t(end),100)), 'g-')
[dEdt,~] = differentiate(Efit,t);

%% Check Amperes law
units = irf_units;
nPlots=5;
for kk=1:nPlots; h(kk) = subplot(nPlots,1,kk); end
isub = 1;

%titlestr = ['n_0 = ' num2str(n0) ',  n_1 = ' num2str(n1) ',  n_2 = ' num2str(n2)];
titlestr = [' '];
w_norm = 1;sqrt(units.mp/units.me);
hca = h(isub); isub=isub+1;
plot(hca,t,mEx,t,Efit(t),'r*',2*pi/w_norm,Efit(2*pi/w_norm),'ko','markerfacecolor','k')
ylabel(hca,'<E_x>')
legend(hca,'E',['E = ' num2str(E0,'%.3f') '*sin(' num2str(omega,'%.3f') '*t)'],'E(t=2*\pi)')
title(hca,titlestr)

if 1
    hca = h(isub); isub=isub+1;
    plot(hca,t,dEdt,t,-4*pi*(mJx0+mJx1+mJx2),'*g')
    legend(hca,'dE/dt','-4\pi*J_{tot}')
    ylabel(hca,'dE/dt')
    title(hca,'Ampere''s law (Gaussian units)')
end

if 0
    hca = h(isub); isub=isub+1;
    plot(hca,t,mJx1)
    ylabel(hca,'<J_{x,1}>')
end

hca = h(isub); isub=isub+1;
plot(hca,t,mJx1,t,mJx2,t,mJx0,t,mJx0+mJx2)
ylabel(hca,'<J_{x}>')
legend(hca,'J_i','J_{e,beam}','J_{e,bg}','J_{e,tot}')

if 0
    hca = h(isub); isub=isub+1;
    plot(hca,t,mu2)
    ylabel(hca,'<v_{x,2}>')
end

if 1
    hca = h(isub); isub=isub+1;
    plot(hca,t,mBx,t,mBy,t,mBz)
    ylabel(hca,'<B>')
    legend(hca,'B_x','B_y','B_z')
end

if 1
    hca = h(isub); isub=isub+1;
    plot(hca,t,mu0*units.me,t,mu1*units.mp,t,mu2*units.me,...
             t,(mu0+mu2)*units.me,t,(mu0+mu2)*units.me+mu1*units.mp,'*',...
             t,(mu0*units.me+mu1*units.mp)*100)
    ylabel(hca,'<p_{x}>')
    legend(hca,'p_{e,bg}','p_i','p_{e,beam}','p_{e,tot}','p_{tot}','p_{e+i,bg}*100')
end






