%% Load PIC
no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');


%% Pick time and location
twpe = 16000;
%no02m.twpelim(twpe).plot_map({'Ey','ne'},'A',1,'sep');
no02m.twpelim(twpe).xlim([97 107]).zlim([-3 3]).plot_map({'Ey','n(4)','ne','Ez','vex'}','A',1,'sep');

%% Plot Ey and contributions
twpe = 17200;
xline_position = no02m.twpelim(twpe).xline_position;
xlim = xline_position(1) + 0.25*[-1 1];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(1*[-1 1]);

varstrs = {{'Ez';'-vex.*By';'vey.*Bx';'divpez';'Ez+vex.*By-vey.*Bx+divpez'},{'Ez';'-vix.*By';'viy.*Bx';'divpiz';'Ez+vix.*By-viy.*Bx-divpiz'},{'Ez';'Jx.*By';'-Jy.*Bx';'divpez';'Jx.*By-Jy.*Bx'}}';
varstrs = {{'Ey';'-vez.*Bx';'vex.*Bz';'-divpey';'Ey+vez.*Bx-vex.*Bz+divpey'},{'vez','-Ey.*Bx./(Bz.^2+Bx.^2)'},{'ne','n(2)','n(4)','n(6)'},{'ni','n(1)','n(3)','n(5)'}}'; % {'vez'},{'Ey.*Bx./(Bz.^2+Bx.^2)'},
varstrs = {{'Ey';'-vez.*Bx';'vex.*Bz';'-divpey';'Ey+vez.*Bx-vex.*Bz+divpey'},{'vez','-Ey.*Bx./(Bz.^2+Bx.^2)'}}';%,{'ne','n(2)','n(4)','n(6)'},{'ni','n(1)','n(3)','n(5)'}}'; % {'vez'},{'Ey.*Bx./(Bz.^2+Bx.^2)'},

varstrs= {{'Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'},{'vez','-Ey.*Bx./(Bz.^2+Bx.^2)'},{'vey','vez','vy(2)','vy([4 6])'},{'jy(2)','jy(4)','jy(6)'},{'n(2)','n(4)','n(6)'}}';
varstrs= {{'Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'},{'vez','-Ey.*Bx./(Bz.^2+Bx.^2)'},{'vey','vez','vy(2)','vy([4 6])'}}';
varstrs= {{'Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'},{'vez','-Ey.*Bx./(Bz.^2+Bx.^2)'},{'vey','Ez.*Bx./(Babs.^2)','(divpiz./ni).*Bx./(Babs.^2)'}}';

h = pic.plot_line('z',varstrs,'smooth',7);
h(2).YLim = [-1 1]*0.99;
c_eval('h(?).FontSize = 15;',1:numel(h))
%%


%%