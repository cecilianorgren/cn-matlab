% See integrate_trajectories
no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');
%trp =
%PICTraj('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/trajectories_paul.h5'); % not yet there

%% Make overview plot to pick starting locations
pic = no02m.twpelim(15000).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);
%pic.plot_map({'vex','vez','vtperp(4)'}','clim',{[-5 5],[-1 1],[0 5]},'A',0.5)
%pic.plot_map({'vex','vez','vtperp(4)','Ez','vey','vy([4 6])','ne','Ez.*Bx./Bx./Bx'}','clim',{[-5 5],[-1 1],[0 5],[-1 1],[-5 5],[-5 5],[0 0.3],[-5 5]},'A',0.5)
%pic.plot_map({'Ez','ne','vey','Ez.*Bx./Bx./Bx','pxy([4 6])','divpy([4 6])','Ey'}','clim',{[-1 1],[0 0.3],[-5 5],[-5 5],0.003*[-1 1],0.3*[-1 1],0.3*[-1 1]},'A',0.5)

%% Reconnection electric field
pic = no02m.twpelim(10000).xlim(mean(no02m.xi)+3*[-5 5]).zlim(3*[-1 1]);
pic = ons.twpelim([13000]).xlim([5 25]).zlim(1.5*[-1 1]);
%varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','n([4 6])','pxy([4 6])','(Bz.*vex-Bx.*vez)','-divpy([4 6])./n([4 6])','Ey','(Bz.*vex-Bx.*vez)-divpy([4 6])./n([4 6])'}';
varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','ne','pexy','(Bz.*vex-Bx.*vez)';'-divpey./ne','Ey','(Bz.*vex-Bx.*vez)-divpey./ne','-vdvey/100','-dveydt/100'}';
clims = {[-5 5],0.5*[-1 1],[0 0.2],0.003*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','Ey','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
clims = {[-5 5],0.5*[-1 1],[0 0.2],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'vey','pexy','Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
clims = {[-5 5],0.01*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

%varstrs = {'n([4 6])./n([2 4 6])','vey','pexy','Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
%clims = {[0 1],7*[-1 1],0.01*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

%varstrs = {'rcz([2])','rcz([4 6])','n([4 6])./n([2 4 6])','vey','Ez'}';
%clims = {[0 0.1],[0 0.1],[0 1],3*[-1 1],[-1 1]};

varstrs = {'vey','Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-divpey./ne-vdvey/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
clims = {[-5 5],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-divpey./ne-vdvey/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
clims = {0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'Ey','(Bz.*vex-Bx.*vez)','pexy','peyz','-divpey./ne','-vdvey/100','-divpey./ne-vdvey/100','Ez'}';
clims = {0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};
clims = {0.1*[-1 1],0.1*[-1 1],0.02*[-1 1],0.02*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],1*[-1 1]};

varstrs = {'Ey','ne','vey','pexy','peyz','teyz','-divpey./ne','-vdvey/100'}';
clims = {0.1*[-1 1],[0 1.2],3.99*[-1 1],0.02*[-1 1],0.02*[-1 1],0.02*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],1*[-1 1]};


varstrs = {'A','Ey','ne','pexy','peyz','-divpey./ne'}';
clims = {[15.7 16.8],0.1*[-1 1],[0 1.2],0.02*[-1 1],0.02*[-1 1],0.1*[-1 1]};



varstrs = {'A','Ey','Bz','vey','Ey+vez.*Bx-vex.*Bz','vex'}';
clims = {[16.2 16.8],0.1*[-1 1],0.1*[-1 1],[-1.5 0],0.02*[-1 1],0.5*[-1 1]};

cmapbr = pic_colors('blue_red');
cmapca = pic_colors('candy4');
cmappa = pic_colors('pasteljet');
cmaps = {cmappa,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};
%varstrs = {'Ez','divpiz./ni','-viy.*Bx+vix.*By'}';
%clims = {1*[-1 1],1*[-1 1],1*[-1 1],[0 0.2]};


varstrs = {'Ey','Bz','vey','vex','-vey.*Ey','Ex','-vey.*Bz'}';
clims = {0.1*[-1 1],0.1*[-1 1],[-1.5 0],0.5*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.5*[-1 1]};
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

h = ons.twpelim([10000 13000]).xlim([5 25]).zlim(0.1*[-1 1]).plot_timemap('xt',varstrs,'clim',clims,'cmap',cmaps,'A',0.01,'smooth',0);
%h = ons.twpelim([10000 13000]).xlim(15+[-0.1 0.1]).zlim(2*[-1 1]).plot_timemap('tz',varstrs,'clim',clims,'cmap',cmaps,'A',0.05);

%h = pic.plot_map(varstrs,'clim',clims,'A',0.1,'sep','smooth',2);
%pic.plot_map(varstrs,'A',0.5,'sep','smooth',2)
%h = no02m.twpelim(10000:1000:20000).xlim(mean(no02m.xi)+3*[-5 5]).zlim(3*[-1 1]).movie(varstrs,'clim',clims,'A',0.5,'sep','smooth',10,'filename',[printpath 'no02m_Ez_balance_ions_smooth2']);
%colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 15;',1:numel(h))

%h = ons.twpelim([10000 13000]).xlim([5 25]).zlim(1.5*[-1 1]).movie(varstrs,'clim',clims,'A',0.1,'sep','smooth',2,'filename',[printpath 'ons_e_eom']);
%ons.twpelim([10000 13000]).xlim([5 25]).zlim(1.5*[-1 1]).movie(varstrs,'clim',clims,'cmap',cmaps,'A',0.1,'sep','smooth',2,'filename',[printpath 'ons_Ay_Ey_Bz']);


%% Nice movies

varstrs = {'vey','vex','Bz','Ey'}';
clims = {2.999*[-1 1],1.999*[-1 1],0.2*[-1 1],0.199*[-1 1]};


cmapbr = pic_colors('blue_red');
cmapca = pic_colors('candy4');
cmappa = pic_colors('pasteljet');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};
ons.twpelim(4000:200:13000).xlim([2 30]).zlim(2*[-1 1]).movie(varstrs,'clim',clims,'cmap',cmaps,'A',0.2,'sep','smooth',2,'filename',[printpath 'no02m_signatures_4']);

%%

varstrs = {'Ey','Bz','vey','vex'}';
clims = {0.1*[-1 1],0.1*[-1 1],[-1.5 0],0.5*[-1 1]};
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};
ons.twpelim(4000:200:13000).xlim([5 25]).zlim(2*[-1 1]).movie(varstrs,'clim',clims,'cmap',cmaps,'A',0.2,'sep','smooth',2,'filename',[printpath 'no02m_signatures_4_5 ']);

%% Background
pic = no02m.twpelim(22000).xlim(mean(no02m.xi)+[-30 30]).zlim(8*[-1 1]);

%h = pic.plot_map({'Ey+(-vex.*Bz+vez.*Bx)'},'clim',{0.3*[-1 1]},'sep','smooth',10);

h = pic.plot_map({'te'},'clim',{[0 0.5]},'sep','smooth',0);
colormap(flipdim(pic_colors('thermal'),1))
c_eval('h(?).FontSize = 15;',1:numel(h))

%% Reconnection signatures
pic = no02m.twpelim(15000).xlim(mean(no02m.xi)+20*[-1 1]).zlim(14*[-0.2 1]);
%varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','n([4 6])','pxy([4 6])','(Bz.*vex-Bx.*vez)','-divpy([4 6])./n([4 6])','Ey','(Bz.*vex-Bx.*vez)-divpy([4 6])./n([4 6])'}';
varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','ne','pexy','(Bz.*vex-Bx.*vez)';'-divpey./ne','Ey','(Bz.*vex-Bx.*vez)-divpey./ne','-vdvey/100','-dveydt/100'}';
clims = {[-5 5],0.5*[-1 1],[0 0.2],0.003*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'vey','(Ez.*Bx-Ex.*Bz)./Babs','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','Ey','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
clims = {[-5 5],0.5*[-1 1],[0 0.2],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'vix','vex','vey','Ez','Jy','Jy.*Ey'}';
clims = {[-5 5],0.01*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

%varstrs = {'n([4 6])./n([2 4 6])','vey','pexy','Ey','(Bz.*vex-Bx.*vez)','-divpey./ne','-vdvey/100','-dveydt/100','(Bz.*vex-Bx.*vez)-divpey./ne-vdvey/100'}';
%clims = {[0 1],7*[-1 1],0.01*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

%varstrs = {'rcz([2])','rcz([4 6])','n([4 6])./n([2 4 6])','vey','Ez'}';
%clims = {[0 0.1],[0 0.1],[0 1],3*[-1 1],[-1 1]};

varstrs = {'vey','vix','vex'}';
clims = {5*[-1 1],0.7*[-1 1],5*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

varstrs = {'Jy','Ey','Jy.*Ey','Ey+vez.*Bx-vex.*Bz','Jy.*(Ey+vez.*Bx-vex.*Bz)'}';
clims = {2*[-1 1],0.3*[-1 1],0.2*[-1 1],0.4*[-1 1],0.3*[-1 1]};


varstrs = {'Jy','Ey','Jy.*Ey','Ey+viz.*Bx-vix.*Bz','Jy.*(Ey+viz.*Bx-vix.*Bz)'}';
clims = {2*[-1 1],0.3*[-1 1],0.2*[-1 1],0.4*[-1 1],0.3*[-1 1]};


%varstrs = {'Jy','Ey','-(viz.*Bx-vix.*Bz)','Ey+viz.*Bx-vix.*Bz'}';
%clims = {2*[-1 1],0.3*[-1 1],0.3*[-1 1],0.4*[-1 1],0.3*[-1 1]};

varstrs = {'Ey','-(viz.*Bx-vix.*Bz)','Ey+viz.*Bx-vix.*Bz'}';
clims = {0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

h = pic.plot_map(varstrs,'clim',clims,'A',0.5,'sep','smooth',2);
%h = pic.plot_map(varstrs,'A',0.5,'sep','smooth',1);
colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 15;',1:numel(h))