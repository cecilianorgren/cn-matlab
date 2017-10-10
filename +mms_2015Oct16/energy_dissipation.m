% Re = E + Ve x B;
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
tintZoom = irf.tint('2015-10-16T10:33:28.00Z/2015-10-16T10:33:31.00Z');
ic = 4;
figure('Name','Energy dissipation')
%c_eval('h = irf_plot({mvaB?,mvaJ?,mvaJe?,mvaJi?,mvaE?.tlim(tintZoom),mvaEVexB?,ne?,RedJ?,RedJe?,RedJi?,EdJ?,EdJe?,EdJi?});',ic)
%c_eval('h = irf_plot({mvaB?,mvaJ?,mvaJe?,mvaJi?,mvaE?.tlim(tintZoom),mvaEVexB?,ne?,RedJ?vec,RedJe?vec,RedJi?vec,EdJ?vec,EdJe?vec,EdJi?vec});',ic)
csys = 'LMN';
switch csys
  case 'LMN'
    c_eval('h = irf_plot({mvaB?,mvaJ?,mvaJe?,mvaJi?,mvaE?.tlim(tintZoom),mvaEVexB?,ne?,mvaRedJ?vec,mvaRedJe?vec,mvaRedJi?vec,mvaEdJ?vec,mvaEdJe?vec,mvaEdJi?vec});',ic)
  case 'GSE'
    c_eval('h = irf_plot({gseB?,gseJ?,gseJe?,gseJi?,gseE?.tlim(tintZoom),gseEVexB?,ne?,RedJ?vec,RedJe?vec,RedJi?vec,EdJ?vec,EdJe?vec,EdJi?vec});',ic)
end

isub = 1;
h(isub).Title.String = [irf_ssub('MMS ?',ic) ' (' csys ')'];
h(isub).YLabel.String = {'B','(nT)'}; isub = isub + 1;
h(isub).YLabel.String = {'J','(nA/m^2)'}; isub = isub + 1;
h(isub).YLabel.String = {'J_e','(nA/m^2)'}; isub = isub + 1;
h(isub).YLabel.String = {'J_i','(nA/m^2)'}; isub = isub + 1;
h(isub).YLabel.String = {'E','(mV/m)'}; isub = isub + 1;
h(isub).YLabel.String = {'E''','(mV/m)'}; isub = isub + 1;
h(isub).YLabel.String = {'n_e','(cm^{-3})'}; isub = isub + 1;
h(isub).YLabel.String = {'E''\cdot J','(nW/m^3)'}; isub = isub + 1;
h(isub).YLabel.String = {'E''\cdot J_e','(nW/m^3)'}; isub = isub + 1;
h(isub).YLabel.String = {'E''\cdot J_i','(nW/m^3)'}; isub = isub + 1;
h(isub).YLabel.String = {'E\cdot J','(nW/m^3)'}; isub = isub + 1;
h(isub).YLabel.String = {'E\cdot J_e','(nW/m^3)'}; isub = isub + 1;
h(isub).YLabel.String = {'E\cdot J_i','(nW/m^3)'}; isub = isub + 1;

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')

for ii = 1:numel(h);
  h(ii).YLabel.Interpreter = 'tex';
  h(ii).FontSize = 10;
end

irf_plot_axis_align

switch csys
  case 'LMN'
    irf_legend(h(1),{'L','M','N'},[0.98 0.95]);
  case 'GSE'
    irf_legend(h(1),{'x','y','z'},[0.98 0.95]);
end

%%