% E/B/phi lower hybrid analysis
tintZoom = irf.tint('2016-12-25T16:00:50.00Z/2016-12-25T16:00:58.00Z');
tintZoomThin = irf.tint('2016-12-25T16:00:55.25Z/2016-12-25T16:00:55.40Z');
tintUTC = tintZoomThin.utc;

[phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(tintZoomThin,gseE1,gseB1scm,gseB1,ne1,'plot',1,'lhfilt',20);

%%
tintZoom = tintZoomThin + 0.05*[-1 1];
x = dirbest;
z = mean(gseB1.tlim(tintZoomThin).data,1); z = z/norm(z);
y = cross(z,x);

vfactor = 3;
dcvPot1_ = irf.ts_scalar(dcvPot1.time,dcvPot1.data(:,1:4));

lmnE1 = gseE1*[x;y;z]'; lmnE1.name = 'fac E';

h=irf_plot({gseB1.tlim(tintZoom),...
            gseB1scm.tlim(tintZoom),...
            gseE1perp.tlim(tintZoom),...
            gseE1par.tlim(tintZoom),...
            lmnE1.tlim(tintZoom),...
            phiEB.*irf.ts_scalar(phiEB.time,repmat([vfactor 1],phiEB.length,1)),...
            -gseVexB1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom)+-dirbest*vbest*vfactor,...
            gseVe1.tlim(tintZoom)-gseVe1.resample(irf_time('2015-11-30T00:24:26.80Z','utc>epochtt')).data,...
            scPot1.tlim(tintZoom),...
            dcvPot1.tlim(tintZoom),...
            ne1.tlim(tintZoom)});
h(1).YLabel.String ='B';          
h(2).YLabel.String ='scm B';          
irf_legend(h(5),{'Ek','E norm','Epar'},[0.95 0.05])
h(6).YLabel.String ='Phi'; irf_legend(h(6),{'phi E','phi B'},[0.95 0.05])
h(7).YLabel.String ='-vexB';
h(8).YLabel.String ='ve';
h(9).YLabel.String ='ve-vph';
h(10).YLabel.String ='ve-veref';
h(11).YLabel.String ='scPot';

add_length_on_top(h(1),vbest*vfactor,1)