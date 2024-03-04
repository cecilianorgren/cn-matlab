p1 = PIC('/Users/cecilia/Data/PIC/varying_tite/data_F1_E05/fields_safety_copy.h5');
p1new = PIC('/Users/cecilia/Data/PIC/varying_tite/data_F1_E05/fields_newdata.h5');
p2 = PIC('/Users/cecilia/Data/PIC/varying_tite/tite_05/fields.h5');
p3 = PIC('/Users/cecilia/Data/PIC/varying_tite/tite_10/fields.h5');

%%
varstrs = {'Ey','Bz','Ez','By','vex','vez','JxBz'}';
clims = {[-1 1],[-1 1],[-2 2],0.5*[-1 1],[-2 2],0.3*[-1 1],0.2*[-1 1]};
cbr = pic_colors('blue_red');
cwa = pic_colors('waterfall');
cmaps = {cbr,cbr,cbr,cbr,cbr,cbr,cbr,cwa};
%p1(50).xlim([10 18]).zlim([0 2]).plot_map(varstrs,'A',0.2,'smooth',2,'clim',clims,'cmap',cmaps,'quiv',{'vex','vez'})
p1(49).xlim([10 18]).zlim([0 2]).plot_map(varstrs,'A',0.2,'smooth',2,'clim',clims,'cmap',cmaps)
%%
varstrs = {'Babs','Ez','By','vex','vez','JxBz./ne','-Jy.*Bx./ne','Jx.*By./ne'}';
clims = {[0 1],[-2 2],0.5*[-1 1],[-2 2],0.3*[-1 1],0.7*[-1 1],0.7*[-1 1],0.7*[-1 1]};
cbr = pic_colors('blue_red');
cwa = pic_colors('waterfall');
cmaps = {cbr,cbr,cbr,cbr,cbr,cbr,cbr,cbr,cbr,cwa};
%p1(50).xlim([10 18]).zlim([0 2]).plot_map(varstrs,'A',0.2,'smooth',2,'clim',clims,'cmap',cmaps,'quiv',{'vex','vez'})
%p1new(p1new.nt).xlim([0 18]).zlim([0 3]).plot_map(varstrs,'A',0.2,'smooth',2,'clim',clims,'cmap',cmaps)
p1new(p1new.nt).plot_map(varstrs,'A',1,'smooth',2,'clim',clims,'cmap',cmaps)
%%
varstrs = {'Ez','By','vperpx([2 4])','vex','vez','vExBx','log10(tepar./teperp)'}';
varstrs = {'Ez','By','vperpx([2 4])','vex','vez','vExBx','JxBz./ne'}';
clims = {[-2 2],0.5*[-1 1],1*[-1 1],2*[-1 1],0.3*[-1 1],1*[-1 1],2*[-1 1]};
cbr = pic_colors('blue_red');
cwa = pic_colors('waterfall');
cmaps = {cbr,cbr,cbr,cbr,cbr,cbr,cbr,cwa,cbr,cbr,cwa};
pic = p1new.twpelim([15000]).xlim([0 25]).zlim([-0.5 2.99]);
%pic.plot_map(varstrs,'A',0.2,'smooth',3,'clim',clims,'cmap',cmaps,'quiv',{'vex','vez'},'quivstep',10)
%pic.movie(varstrs,'A',0.1,'smooth',3,'clim',clims,'cmap',cmaps,'filename',[printpath 'rerec3'])
pic.movie(varstrs,'A',0.1,'smooth',3,'clim',clims,'cmap',cmaps,'filename',[printpath 'rerec7_newdata'])

%% plotline, Ohm's law
comp = 'x';
twpe = [1000:1000:9000];
twpe = 9450;
xlim = [11 18];
zlim = 0.385+0.05*[-1 1];
%pic = no02.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
pic = p1.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'ni'};{'divpx([1 3 5])','divpx([1])','divpx([3 5])','divpx(1+[1 3 5])'};{'dvxdt([1])','vdvx(1)'};{'vex','vix','vExBx'};{'t([1 3 5])'}};
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'divpx([1 3 5])','vxBx([1 3 5])','Ex','-divpx([1 3 5])+vxBx([1 3 5])'};{'dvxdt([1])','vdvx(1)'}};
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'divpx([1 3 5])','vxBx([1 3 5])','Ex','-divpx([1 3 5])+vxBx([1 3 5])'};{'t([1 3 5])'}};
varstrs = {{'Bz','Ey','Jy'};   {'pi','pe'};{'JxBx','JxBy','JxBz'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','-vixBx+divpix'};{'Ex','-vexBx','-divpex','-vexBx-divpex'};{'Ex+vixBx','Ex+vxBx(1)','Ex+vxBx([3 5])','Ex+vexBx'}};
varstrs = {{'Bz','Ey','Jy','pi','pe'};{'jiy','jey'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','-vixBx+divpix'};{'Ex','-vexBx','-divpex','-vexBx-divpex'};{'Ex+vixBx','Ex+vexBx'}};
%varstrs = {{'Bz','Ey','Jy','pi','pe'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','dvixdt','vdvix','-vixBx+divpix+dvixdt+vdvix'};{'Ex','-vexBx','-divpex','dvexdt','vdvex','-vexBx-divpex'}};     
%varstrs = {{'Bz','Ey','Jy'};{'Jx','Jy','Jz'};{'JxBx','JxBy','JxBz'};{'Ey','vixBy','divpey','JxBy','JxBy-vixBy+divpey'};{'Ey','divpiy','vixBy','-vixBy+divpiy'};{'Ey','-vexBy','-divpey','-vexBy-divpey'};{'Ey+vixBy','Ey+vexBy'}};
%varstrs = {{'Bz','Ey','Jy'};   {'JxBx','JxBy','JxBz'};{'Ey','vixBy','divpey','JxBy','JxBy-vixBy+divpey'};{'Ey','divpiy','vixBy','-vixBy+divpiy'};{'Ey','-vexBy','-divpey','-vexBy-divpey'};{'Ey+vixBy','Ey+vexBy'}};
%varstrs = {{'Bz','Ey'};{'Jz','Jy','Jz'};{'JxBx','JxBy','JxBz'};{'Ex','JxBx','vxBx([1 3])','divpx(1+[1 3])','JxBx-vxBx([1 3])+divpx(1+[1 3])'};{'Ex','divpx([1 3])','vxBx([1 3])','-vxBx([1 3])+divpx([1 3])'};{'Ex','-vxBx(1+[1 3])','-divpx(1+[1 3])','-vxBx(1+[1 3])-divpx(1+[1 3])'}};
%varstrs = {{'Bz'};{'t([1 3 5])'}};
varstrs = {{'By','Ez','Jx'};...
  {'jix','jex'};...
  {'(Jx.*By-Jy.*Bx)./ni','Jx.*By./ni','-Jy.*Bx./ni'};
  {'Ez','vix.*By','-viy.*Bx','Ez+vix.*By-viy.*Bx'};...
  {'Ez','vex.*By','-vey.*Bx','Ez+vex.*By-vey.*Bx'};...
  {'Ez','vix.*By-viy.*Bx','(Jx.*By-Jy.*Bx)./ni','-divpex./ne'};...
  {'Ez+vix.*By-viy.*Bx','(Jx.*By-Jy.*Bx)./ni','-divpex./ne'};...
  {'Ez+vix.*By-viy.*Bx','-divpix./ni'};...
  {'Ez+vex.*By-vey.*Bx','divpex./ne'}};%;...
  %{'Ez','vix.*By','viy.*Bx','-divpex./ne','Jx.*By./ni','-Jy.*Bx./ni'}};%;...
  %{'Ez','vixBx','divpix','-vixBx+divpix'};...
  %{'Ez','-vexBx','-divpex','-vexBx-divpex'};...
  %{'Ez+vixBx','Ex+vexBx'}};
h = pic.plot_line(comp,varstrs,'smooth',2);

%%

zmax = 12.5;
l = 2;
B0 = 1;
flux0 = B0*l*log(cosh(zmax/l));
rel_flux_add = 1;
flux_add = flux0*rel_flux_add;
Epeak = 0.5;

E0 = Epeak*exp(1);
ts = flux_add/E0;
Ed = @(t,ts,E0) E0*(t/ts).*exp(-t/ts);
Eypeak = Ed(ts,ts,E0);


t = p1.twci;
t_ = linspace(0,t(end),100);
Ey_sim = squeeze(mean(mean(p1.xlim(mean(p1.xi)+[-0.2 0.2]).zlim(max(p1.zi)+[-0.2 0]).Ey)));
intE_sim = cumtrapz(t,Ey_sim);
intE_ana = cumtrapz(t_,Ed(t_,ts,E0));
hca = subplot(2,1,1);
plot(hca,t,Ey_sim,...
     t_,Ed(t_,ts,E0),'--',...
     t,Ed(t,ts,E0),...
     ts*[1 1],[0 Epeak],'k--');
   
hca = subplot(2,1,2);
plot(hca,...
     t,intE_sim,...
     t_,intE_ana,...
     t_,t_*0+flux_add,'k-',...
     t_,t_*0+flux0,'r--',...
     ts*[1 1],[0 flux_add],'k--');