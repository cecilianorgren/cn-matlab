% use times from sql database to fetch and work with mms data from spis.

%% Get times
% eventId is start of burst interval in format YYYYMMDD_hhmmss
query = ['select eventId,X,Y,Z,Bmin,t_0,t_12,t_88,Lx,Ly,Lz,Mx,My,Mz,Nx,Ny,Nz from MMS3'];
[eventId,X,Y,Z,Bmin,t_0,t_12,t_88,Lx,Ly,Lz,Mx,My,Mz,Nx,Ny,Nz] = mysql(query); %SymH assumed equal to DST
clear query;

query = ['select eventId,Vex_tmaxJ,Vey_tmaxJ,Vez_tmaxJ from MMS3_elctrns'];
[eventId_elctrns,Vex_tmaxJ,Vey_tmaxJ,Vez_tmaxJ] = mysql(query); %SymH assumed equal to DST
% get lmn velocities

[C,ia,ib] = intersect(eventId,eventId_elctrns);
VeL_tmaxJ = Vex_tmaxJ(ib).*Lx(ia) + Vey_tmaxJ(ib).*Ly(ia) + Vez_tmaxJ(ib).*Lz(ia);
VeM_tmaxJ = Vex_tmaxJ(ib).*Mx(ia) + Vey_tmaxJ(ib).*My(ia) + Vez_tmaxJ(ib).*Mz(ia);
VeN_tmaxJ = Vex_tmaxJ(ib).*Nx(ia) + Vey_tmaxJ(ib).*Ny(ia) + Vez_tmaxJ(ib).*Nz(ia);

%cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_0,'UniformOutput', false)

% Make EpochTT list
t_0_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_0,'UniformOutput', false)));
t_12_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_12,'UniformOutput', false)));
t_88_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_88,'UniformOutput', false)));

dt = t_88_epochtt-t_12_epochtt;

if 1
  %%
  hca = subplot(1,1,1);%h(isub); isub = isub + 1;
  edgesVeL_tmaxJ = linspace(min(VeL_tmaxJ),max(VeL_tmaxJ),50);
  edgesVeM_tmaxJ = linspace(min(VeM_tmaxJ),max(VeM_tmaxJ),50);
  edgesVeN_tmaxJ = linspace(min(VeN_tmaxJ),max(VeN_tmaxJ),50);
  
  [count edges mid loc] = histcn([VeL_tmaxJ,VeM_tmaxJ], edgesVeL_tmaxJ, edgesVeM_tmaxJ);
  pcolor(hca,mid{1:2},count');
  hca.XLabel.String = 'VeM_ tmaxJ';
  hca.YLabel.String = 'VeM_ tmaxJ';
end
