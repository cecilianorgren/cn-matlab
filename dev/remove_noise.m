nMean = [5,3,3,3]; nThresh = 3;

c_eval('PD_counts = iPDist?_counts;',ic)
c_eval('PD_orig = iPDist?;',ic)
PD_clean = PD_orig.remove_noise(nMean,nThresh,PD_counts);
PD_diff = PD_orig+-PD_clean;

counts = 5;
nMovMean = [7 7];
matEmask = PD_counts.find_low_counts('counts',counts,'nMovMean',nMovMean,'output','mat');
%emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit
tsElow = PD_counts.find_low_counts('counts',counts,'nMovMean',nMovMean,'output','energy');

PD_orig_notmasked = PD_orig;
PD_clean_notmasked = PD_clean;
PD_diff_notmasked = PD_diff;

PD_orig_mask = PD_orig.mask('energy','mat',matEmask);
PD_clean_mask = PD_clean.mask('energy','mat',matEmask);
PD_diff_mask = PD_diff.mask('energy','mat',matEmask);

h = irf_plot({PD_orig.deflux.omni.specrec,PD_clean_notmasked.deflux.omni.specrec,PD_clean_mask.deflux.omni.specrec}); c_eval('h(?).YScale = ''log'';',1:numel(h))
hlinks = linkprop(h,{'CLim','YLim'});


PD_orig = PD_orig_mask;
PD_clean = PD_clean_mask;
PD_diff = PD_diff_mask;