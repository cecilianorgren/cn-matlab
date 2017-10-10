a=rand(100,1)*100;
plot(a)
nanmean(log10(a))
log10(nanmean(a))
%%
cd /Users/Cecilia/Data/BM/20070831
t1=[2007 08 31 10 19 06.90]; t2=[2007 08 31 10 19 07.50]; % lhdw
tint=toepoch([t1;t2])';
%data_structure = c_caa_distribution_data('C3_CP_PEA_3DXPH_PSD');
for k=1:8
    h(k)=subplot(2,4,k);
end
c_caa_plot_distribution_function(h(1),'tint',tint(1),'t_display','tags','polar',data_structure,data_structure);
c_caa_plot_distribution_function(h(2),'tint',tint(1),'t_display','given','polar',data_structure,data_structure);
c_caa_plot_distribution_function(h(3),'tint',tint(1),'t_display','tags','polar',data_structure);
c_caa_plot_distribution_function(h(4),'tint',tint(1),'t_display','given','polar',data_structure);
c_caa_plot_distribution_function(h(5),'tint',tint,'t_display','tags','polar',data_structure,data_structure);
c_caa_plot_distribution_function(h(6),'tint',tint,'t_display','given','polar',data_structure,data_structure);
c_caa_plot_distribution_function(h(7),'tint',tint,'t_display','tags','polar',data_structure);
c_caa_plot_distribution_function(h(8),'tint',tint,'t_display','given','polar',data_structure);

%c_caa_plot_distribution_function(h2,'tint',tint,'t_display','tags','polar',data_structure);
%%
a='1';
b='0';
switch [a,b]
    case '12'
        1
    case '10'
        2
end
%%
