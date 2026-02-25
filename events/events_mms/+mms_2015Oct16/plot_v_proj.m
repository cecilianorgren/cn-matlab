% plot_v_proj

t1 = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z'); t1 = t1(1);

nrows = 3;
ncols = 4;
isub = 1;
figure;
for kk = 1:nrows;
  for pp = 1:ncols;        
    h(isub) = subplot(nrows,ncols,isub); hca = h(isub); 
    tplot = t1 + isub; 
    tints(isub) = tplot;
    fpi_plot_proj(disDist1,tplot,'xz')
    isub = isub + 1;
  end
end

%%
isub = 1;
for kk = 1:nrows*ncols;
  axis(h(kk),'square')
  axis(h(kk),'equal')
  h(kk).YLim = 500*[-1 1];
  h(kk).XLim = 500*[-1 1];
end
