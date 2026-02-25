
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% NB: all input should be in same coordinate system
ePDist1_model = mms.make_model_dist(ePDist1,dmpaB1,scPot1,ne1,dbcsVe1,gseTe1);

ePDist1_remerr = ePDist1;
c_eval('ePDist?_remerr.data(ePDist?_remerr.data<ePDistErr?.data*1.1) = 0;',1)
ierr = ePDist1.data==ePDistErr1.data;
ierr_ = ePDist1.data<ePDistErr1.data*1.1;
%%
h = irf_plot({...
  ePDist1.convertto('s^3/km^6').omni.specrec('en'),...
  ePDistErr1.convertto('s^3/km^6').omni.specrec('en'),...
  ePDist1_remerr.convertto('s^3/km^6').omni.specrec('en'),...
  ePDist1_model.omni.specrec('en')});

for ip = 1:numel(h)
  h(ip).YScale = 'log';
  colorbar('peer',h(ip))  
  h(ip).CLim = [-6 6];
end

irf_plot_axis_align