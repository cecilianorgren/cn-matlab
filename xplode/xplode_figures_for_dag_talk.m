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
pic = no02m.twpelim(14000).xlim(mean(no02m.xi)+[-4 4]).zlim([-1 1]);
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


h = pic.plot_map(varstrs,'clim',clims,'A',0.2,'sep','smooth',0);
%pic.plot_map(varstrs,'A',0.5,'sep','smooth',2)
colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 15;',1:numel(h))