function tred_plot_cores(core,satellite,position)
% TRED_PLOT_CORES
%
%
%

axis([floor(min(position(:,1))) ceil(max(position(:,1))),...
      floor(min(position(:,2))) ceil(max(position(:,2))),...
      floor(min(position(:,3))) ceil(max(position(:,3)))]);
  
str_core_satellite=cellstr([num2str(core),repmat('-',max(size(core)),1),num2str(satellite)]);
text(position(:,1),position(:,2),position(:,3),str_core_satellite);

xlabel('x');ylabel('y');zlabel('z')