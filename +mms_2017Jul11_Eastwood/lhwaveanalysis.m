% E/B/phi lower hybrid analysis
tintZoom = irf.tint('2017-07-11T22:34:01.00Z/2017-07-11T22:34:03.00Z'); %20151112071854
%tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
%tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
%tintZoomThin = irf.tint('2015-11-30T00:24:26.68Z/2015-11-30T00:24:26.75Z');
%tintZoomThin = irf.tint('2015-11-30T00:22:36.50Z/2015-11-30T00:22:37.50Z');
%tintZoomThin = irf.tint('2015-11-30T00:23:16.00Z/2015-11-30T00:23:17.00Z');
tintLH = tintZoom;
tintUTC = tintLH.utc;
ffilt = 3;
mmsid = 1;
c_eval('[phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(tintLH,gseE?,gseB?scm,gseB?,ne?,''plot'',1,''lhfilt'',ffilt);',mmsid)

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