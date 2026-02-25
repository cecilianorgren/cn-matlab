% plot only isolated energy sweeps
% new plot with more pitch angles
tstart = toepoch([2007 08 31 10 17 30.00]);
tstop  = toepoch([2007 08 31 10 18 15.00]);
dt     = 0.125; % one energy sweep interval

if ~exist('peaPSD3','var'); load peaPSD; end
if ~exist('diE3','var'); load matlabdiEB; end

t = tstart;
nPlot = 1;
figure('units','normalized','position',[0 0 1 1])
while t < tstop
    for k = 1:8
    h(k) = subplot(2,4,k); 
    hold(h(k),'on')
    c_caa_plot_distribution_function(...
        h(k),...
        'tint',t,...
        'cross-section',...
        'pitchangle',peaPSD3.f_cs,...
        't_display','tags',...
        peaPSD3);    
    c_caa_plot_distribution_function(...
        h(k),...
        'tint',t,...
        'cross-section',...
        'pitchangle',peaPSD4.f_cs,...
        't_display','tags',...
        peaPSD4);   
    set(h(k),'ylim',[1e-5 1e2],'xlim',[4e1 3e4],'yscale','log','xscale','log')
    %hold(h(k),'off')
    t = t + dt;
    end
    cn.print(['energy_sweeps_' num2str(nPlot)],'EnergySweeps2x4_C3C4')
    nPlot = nPlot + 1;
    delete(h);
end

%tint = [toepoch(t1) toepoch(t2)];

