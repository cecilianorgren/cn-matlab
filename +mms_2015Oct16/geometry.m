% Geometry
tint = irf.tint('2015-10-16T10:32:30.00Z/2015-10-16T10:34:10.00Z');

%% LMN coordinate system from minimum variance analysis of B
[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
L = v(1,:); M = v(2,:); N = v(3,:);

%% Normal to magnetospheric boundary layer
% Normal to magnetospheric boundary: Vector normal to line between mms3 and mms4
gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
gseM = irf.ts_vec_xyz(gseVec34.time,M);
gseNorm34 = gseVec34.cross(gseM);
gseNormalMSP = gseNorm34/gseNorm34.abs;

mspN = -gseNormalMSP.data;
mspM = M;
mspL = cross(mspM,mspN);

mvaVec34 = irf.ts_vec_xyz(gseVec34.time,[gseVec34.dot(L).data gseVec34.dot(M).data gseVec34.dot(N).data]);
mvaM = irf.ts_vec_xyz(mvaVec34.time,[0 1 0]);
mvaNorm34 = mvaVec34.cross(mvaM);
mvaNormalMSP = mvaNorm34/mvaNorm34.abs;

% Normal to magnetosheath boundary: Vector normal to line between mms1 and mms4
gseVec14 = gseR4-gseR1; gseVec14 = gseVec14.resample(tint.start);
gseM = irf.ts_vec_xyz(gseVec14.time,M);
gseNorm14 = gseVec14.cross(gseM);
gseNormalMSH = gseNorm14/gseNorm14.abs;

mshN = -gseNormalMSH.data;
mshM = M;
mshL = cross(mshM,mshN);

mvaVec14 = irf.ts_vec_xyz(gseVec14.time,[gseVec14.dot(L).data gseVec14.dot(M).data gseVec14.dot(N).data]);
mvaM = irf.ts_vec_xyz(mvaVec14.time,[0 1 0]);
mvaNorm14 = mvaVec14.cross(mvaM);
mvaNormalMSH = mvaNorm14/mvaNorm14.abs;

[mspL;mspM;mspN]
[L;M;N]
[mshL;mshM;mshN]

hca = subplot(1,1,1); hold(hca,'on')
plot_quivers(hca,mspL,[0 0 0],'k')
plot_quivers(hca,mspM,[0 0 0],'k')
plot_quivers(hca,mspN,[0 0 0],'k')
plot_quivers(hca,L,[0 0 0],'b')
plot_quivers(hca,M,[0 0 0],'b')
plot_quivers(hca,N,[0 0 0],'b')
plot_quivers(hca,mshL,[0 0 0],'r')
plot_quivers(hca,mshM,[0 0 0],'r')
plot_quivers(hca,mshN,[0 0 0],'r')
axis(hca,'equal')

hold(hca,'on')


%% Compare different quantities in different coordinate systems
%% Set up coordinate system
% N: minimum variance of B
[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
mvaL = v(1,:); mvaM = v(2,:); mvaN = v(3,:);

% N: minimum variance of J
[out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-10-16T10:33:23.000Z/2015-10-16T10:33:32.000Z')));
mvajL = -v(2,:); mvajM = -v(1,:); mvajN = v(3,:);

% N: maximum variance of E    
tint = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.400Z');
[out,l,v] = irf_minvar(gseE3.tlim(tint));
mvaeL = v(2,:); mvaeM = -v(3,:); mvaeN = -v(1,:);

% N: magnetosheath side normal derived from mms1 and mms4
gseVec14 = gseR4-gseR1; gseVec14 = gseVec14.resample(tint.start);
gseM = irf.ts_vec_xyz(gseVec14.time,M);
gseNorm14 = gseVec14.cross(gseM);
gseNormalMSH = gseNorm14/gseNorm14.abs;
mshN = -gseNormalMSH.data;
mshM = M;
mshL = cross(mshM,mshN);

% N: magnetosphere side normal derived from mms1 and mms4
gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
gseM = irf.ts_vec_xyz(gseVec34.time,M);
gseNorm34 = gseVec34.cross(gseM);
gseNormalMSP = gseNorm34/gseNorm34.abs;
mspN = -gseNormalMSP.data;
mspM = M;
mspL = cross(mspM,mspN);
%% Rotate data
gradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4);

c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(mvaL).data gseR?.dot(mvaM).data gseR?.dot(mvaN).data]);')
c_eval('mvaB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(mvaL).data gseB?.dot(mvaM).data gseB?.dot(mvaN).data]);')
c_eval('mvaE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(mvaL).data gseE?.dot(mvaM).data gseE?.dot(mvaN).data]);')
c_eval('mvaVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(mvaL).data gseVe?.dot(mvaM).data gseVe?.dot(mvaN).data]);')
c_eval('mvaVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(mvaL).data gseVi?.dot(mvaM).data gseVi?.dot(mvaN).data]);')
c_eval('mvaJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(mvaL).data gseJ?.dot(mvaM).data gseJ?.dot(mvaN).data]); mvaJ?.units = gseJ?.units;')
c_eval('mvaPe? = mms.rotate_tensor(gsePe?,''rot'',mvaL,mvaM,mvaN); mvaPe? = irf.ts_tensor_xyz(mvaPe?.time,mvaPe?.data); mvaPe?.units = Pe?.units;',ic)
mvaGradPe = irf.ts_vec_xyz(gradPe.time,[gradPe.dot(mvaL).data gradPe.dot(mvaM).data gradPe.dot(mvaN).data]);
c_eval('mvaVexB? = mvaVe?.resample(mvaB?.time).cross(mvaB?)*1e-3; mvaVexB?.units = ''mV/m'';')
c_eval('mvaVixB? = mvaVi?.resample(mvaB?.time).cross(mvaB?)*1e-3; mvaVixB?.units = ''mV/m'';')
c_eval('mvaJxB? = mvaJ?.resample(mvaB?.time).cross(mvaB?);')

c_eval('mspR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(mspL).data gseR?.dot(mspM).data gseR?.dot(mspN).data]);')
c_eval('mspB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(mspL).data gseB?.dot(mspM).data gseB?.dot(mspN).data]);')
c_eval('mspE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(mspL).data gseE?.dot(mspM).data gseE?.dot(mspN).data]);')
c_eval('mspVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(mspL).data gseVe?.dot(mspM).data gseVe?.dot(mspN).data]);')
c_eval('mspVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(mspL).data gseVi?.dot(mspM).data gseVi?.dot(mspN).data]);')
c_eval('mspJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(mspL).data gseJ?.dot(mspM).data gseJ?.dot(mspN).data]); mspJ?.units = gseJ?.units;')
c_eval('mspPe? = mms.rotate_tensor(gsePe?,''rot'',mspL,mspM,mspN); mspPe? = irf.ts_tensor_xyz(mspPe?.time,mspPe?.data); mspPe?.units = Pe?.units;',ic)
mspGradPe = irf.ts_vec_xyz(gradPe.time,[gradPe.dot(mspL).data gradPe.dot(mspM).data gradPe.dot(mspN).data]);
c_eval('mspVexB? = mspVe?.resample(mspB?.time).cross(mspB?)*1e-3; mspVexB?.units = ''mV/m'';')
c_eval('mspVixB? = mspVi?.resample(mspB?.time).cross(mspB?)*1e-3; mspVixB?.units = ''mV/m'';')
c_eval('mspJxB? = mspJ?.resample(mspB?.time).cross(mspB?);')

c_eval('mshR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(mshL).data gseR?.dot(mshM).data gseR?.dot(mshN).data]);')
c_eval('mshB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(mshL).data gseB?.dot(mshM).data gseB?.dot(mshN).data]);')
c_eval('mshE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(mshL).data gseE?.dot(mshM).data gseE?.dot(mshN).data]);')
c_eval('mshVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(mshL).data gseVe?.dot(mshM).data gseVe?.dot(mshN).data]);')
c_eval('mshVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(mshL).data gseVi?.dot(mshM).data gseVi?.dot(mshN).data]);')
c_eval('mshJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(mshL).data gseJ?.dot(mshM).data gseJ?.dot(mshN).data]); mshJ?.units = gseJ?.units;')
c_eval('mshPe? = mms.rotate_tensor(gsePe?,''rot'',mshL,mshM,mshN); mshPe? = irf.ts_tensor_xyz(mshPe?.time,mshPe?.data); mshPe?.units = Pe?.units;',ic)
mshGradPe = irf.ts_vec_xyz(gradPe.time,[gradPe.dot(mshL).data gradPe.dot(mshM).data gradPe.dot(mshN).data]);
c_eval('mshVexB? = mshVe?.resample(mshB?.time).cross(mshB?)*1e-3; mshVexB?.units = ''mV/m'';')
c_eval('mshVixB? = mshVi?.resample(mshB?.time).cross(mshB?)*1e-3; mshVixB?.units = ''mV/m'';')
c_eval('mshJxB? = mshJ?.resample(mshB?.time).cross(mshB?);')
%% Example 3 (Ohm's law, 4 sc): Initialize plot
[h1,h2] = initialize_combined_plot(10,3,2,3,'horizontal');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');

%% Example 3 (Ohm's law, 4 sc): Plot timeseries with mva of B
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:32.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeM
  hca = irf_panel('Ve M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).y,mvaVe2.tlim(tint).y,mvaVe3.tlim(tint).y,mvaVe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'v_e_M','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeN
  hca = irf_panel('Ve N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).z,mvaVe2.tlim(tint).z,mvaVe3.tlim(tint).z,mvaVe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'v_e_N','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).x,mvaE2.tlim(tint).x,mvaE3.tlim(tint).x,mvaE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).y,mvaE2.tlim(tint).y,mvaE3.tlim(tint).y,mvaE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB N
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.x+1*mvaVexB1.resample(mvaE1.time).x,...
                mvaE2.x+1*mvaVexB2.resample(mvaE2.time).x,...
                mvaE3.x+1*mvaVexB3.resample(mvaE3.time).x,...
                mvaE4.x+1*mvaVexB4.resample(mvaE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 0 % E + vexB N
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.y+1*mvaVexB1.resample(mvaE1.time).y,...
                mvaE2.y+1*mvaVexB2.resample(mvaE2.time).y,...
                mvaE3.y+1*mvaVexB3.resample(mvaE3.time).y,...
                mvaE4.y+1*mvaVexB4.resample(mvaE4.time).y},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 1 % E + vexB N
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align

isub = 1;         
hca = h2(isub); isub = isub + 1; 
lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 -1 0])
axis(hca,'square')

hca = h2(isub); isub = isub + 1; 
plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 0 1])
axis(hca,'square')
%% Example 3 (Ohm's law, 4 sc): Plot timeseries with msp normal from mms3 and mms4
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:32.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspB1.x.tlim(tint),mspB2.x.tlim(tint),mspB3.x.tlim(tint),mspB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspB1.y.tlim(tint),mspB2.y.tlim(tint),mspB3.y.tlim(tint),mspB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspB1.z.tlim(tint),mspB2.z.tlim(tint),mspB3.z.tlim(tint),mspB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspVe1.tlim(tint).x,mspVe2.tlim(tint).x,mspVe3.tlim(tint).x,mspVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeM
  hca = irf_panel('Ve M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspVe1.tlim(tint).y,mspVe2.tlim(tint).y,mspVe3.tlim(tint).y,mspVe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'v_e_M','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeN
  hca = irf_panel('Ve N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspVe1.tlim(tint).z,mspVe2.tlim(tint).z,mspVe3.tlim(tint).z,mspVe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'v_e_N','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspE1.tlim(tint).x,mspE2.tlim(tint).x,mspE3.tlim(tint).x,mspE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspE1.tlim(tint).y,mspE2.tlim(tint).y,mspE3.tlim(tint).y,mspE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mspE1.tlim(tint).z,mspE2.tlim(tint).z,mspE3.tlim(tint).z,mspE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB N
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mspE1.x+1*mspVexB1.resample(mspE1.time).x,...
                mspE2.x+1*mspVexB2.resample(mspE2.time).x,...
                mspE3.x+1*mspVexB3.resample(mspE3.time).x,...
                mspE4.x+1*mspVexB4.resample(mspE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 0 % E + vexB N
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mspE1.y+1*mspVexB1.resample(mspE1.time).y,...
                mspE2.y+1*mspVexB2.resample(mspE2.time).y,...
                mspE3.y+1*mspVexB3.resample(mspE3.time).y,...
                mspE4.y+1*mspVexB4.resample(mspE4.time).y},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 1 % E + vexB N
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mspE1.z+1*mspVexB1.resample(mspE1.time).z,...
                mspE2.z+1*mspVexB2.resample(mspE2.time).z,...
                mspE3.z+1*mspVexB3.resample(mspE3.time).z,...
                mspE4.z+1*mspVexB4.resample(mspE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align

isub = 1;         
hca = h2(isub); isub = isub + 1; 
lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
plot_lmn3D(hca,mspR1,mspR2,mspR3,mspR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 -1 0])
axis(hca,'square')

hca = h2(isub); isub = isub + 1; 
plot_lmn3D(hca,mspR1,mspR2,mspR3,mspR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 0 1])
axis(hca,'square')
%% Example 3 (Ohm's law, 4 sc): Plot timeseries with msh normal from mms1 and mms4
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:32.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshB1.x.tlim(tint),mshB2.x.tlim(tint),mshB3.x.tlim(tint),mshB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshB1.y.tlim(tint),mshB2.y.tlim(tint),mshB3.y.tlim(tint),mshB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshB1.z.tlim(tint),mshB2.z.tlim(tint),mshB3.z.tlim(tint),mshB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshVe1.tlim(tint).x,mshVe2.tlim(tint).x,mshVe3.tlim(tint).x,mshVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeM
  hca = irf_panel('Ve M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshVe1.tlim(tint).y,mshVe2.tlim(tint).y,mshVe3.tlim(tint).y,mshVe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'v_e_M','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % VeN
  hca = irf_panel('Ve N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshVe1.tlim(tint).z,mshVe2.tlim(tint).z,mshVe3.tlim(tint).z,mshVe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'v_e_N','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshE1.tlim(tint).x,mshE2.tlim(tint).x,mshE3.tlim(tint).x,mshE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshE1.tlim(tint).y,mshE2.tlim(tint).y,mshE3.tlim(tint).y,mshE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mshE1.tlim(tint).z,mshE2.tlim(tint).z,mshE3.tlim(tint).z,mshE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB L
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mshE1.x+1*mshVexB1.resample(mshE1.time).x,...
                mshE2.x+1*mshVexB2.resample(mshE2.time).x,...
                mshE3.x+1*mshVexB3.resample(mshE3.time).x,...
                mshE4.x+1*mshVexB4.resample(mshE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 0 % E + vexB M
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mshE1.y+1*mshVexB1.resample(mshE1.time).y,...
                mshE2.y+1*mshVexB2.resample(mshE2.time).y,...
                mshE3.y+1*mshVexB3.resample(mshE3.time).y,...
                mshE4.y+1*mshVexB4.resample(mshE4.time).y},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end
if 1 % E + vexB N
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mshE1.z+1*mshVexB1.resample(mshE1.time).z,...
                mshE2.z+1*mshVexB2.resample(mshE2.time).z,...
                mshE3.z+1*mshVexB3.resample(mshE3.time).z,...
                mshE4.z+1*mshVexB4.resample(mshE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  hca.YLim = [-10 10];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align

isub = 1;         
hca = h2(isub); isub = isub + 1; 
lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
plot_lmn3D(hca,mshR1,mshR2,mshR3,mshR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 -1 0])
axis(hca,'square')

hca = h2(isub); isub = isub + 1; 
plot_lmn3D(hca,mshR1,mshR2,mshR3,mshR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'})
hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
view(hca,[0 0 1])
axis(hca,'square')

%% Example 3 (Ohm's law, 4 sc): Plot 2x sc positions + 4x projection, electrons 

tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:34:00.00Z');
%tt = irf.tint('2015-10-16T10:33:26.350Z',0.1);
%tt = irf.tint('2015-10-16T10:33:30.440Z',0.1)+(-0.18);
%tt = irf.tint('2015-10-16T10:33:27.920Z',0.1);
%tt=tt(1);

times = desDist1.tlim(tintZoom).time;


for it = 119;1:times.length  % 644 - ring distributions
  time = times(it);
  tint = time;
  tt = time;
  
  
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,tt.epochUnix','green');

  vectors = {'B',{mvaB1.resample(tt(1)),mvaB2.resample(tt(1)),mvaB3.resample(tt(1)),mvaB4.resample(tt(1))},2;...
             'Ve',{mvaVe1.resample(tt(1)),mvaVe2.resample(tt(1)),mvaVe3.resample(tt(1)),mvaVe4.resample(tt(1))},70};        

  isub = 1;         
  hca = h2(isub); isub = isub + 1; 
  lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 -1 0])
  axis(hca,'square')

  hca = h2(isub); isub = isub + 1; 
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 0 1])
  axis(hca,'square')

  %for ic = 1:4;
  %hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,desDist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  vlim = 15*1e3;
    elevlim = 20;
    strCMap = 'jet';
    %energies =  [30 220];
    projclim = [0 4.5];
    palim = [1e-3 1e6];
    skymapEnergy = [65 278];


  for ic = 1:4;
      c_eval('dist = desDist?;',ic)

      c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
      hatVe0 = double(irf_norm(Ve0));    

      % Get mean magnetic field direction
      c_eval('B0 = gseB?.resample(time).data;',ic); 
      hatB0 = double(irf_norm(B0));
      c_eval('E0 = gseE?.resample(time).data;',ic); 
      hatE0 = double(irf_norm(E0));
      hatExB0 = cross(hatE0,hatB0);

      vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

      % Projection coordinate system
      if 1
        x = hatB0;
        y = hatExB0;
        z = cross(x,y);
      else
        x = [1 0 0];
        y = [0 1 0];
        z = [0 0 1];
      end



      %isub = ic;

      if 0 % Plot psd 0 90 180
        hca = h2(isub); isub = isub + 1;
        psdtint = tint;%+[-1 1]*0.03;
        c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
        TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
        %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
        %hca.Title.String = irf_ssub('C?',ic);
        hca.Title.String = TitleStr;
      end

      if 0 % Plot skymap for a given energy
        hca = h2(isub); isub = isub + 1;      
        c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
        %hca.Title.String = hca.Title.String{2};

        % Plot skymap for a given energy
        hca = h2(isub); isub = isub + 1;      
        c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
        %hca.Title.String = hca.Title.String{2};
      end

      % Plot project ion onto a plane
      if 0
        hca = h2(isub); isub = isub + 1; 
        mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
        hca.Title.String = '';
        colormap(hca,strCMap)
      end

      hca = h2(isub); isub = isub + 1; 
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
      hca.Title.String = titleStr;
      colormap(hca,strCMap)

      if 0
        hca = h2(isub); isub = isub + 1; 
        mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
        hca.Title.String = '';
        colormap(hca,strCMap)
      end
    end
    pause(1)
    %cn.print([irf_ssub('ohm_sc_pos_4sc',1) irf_time(time,'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end
  
if 0
hca = h2(3);
tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
c_eval('tts? = double((mvaB?.tlim(tintQuiver).time.ttns-mvaB?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
c_eval('R? = mvaR?-mvaR0;')
velocity = [0 0 30];
c_eval('xyz? = irf.ts_vec_xyz(mvaB?.tlim(tintQuiver).time,tts?*velocity);')
c_eval('quivers? = mvaB?.tlim(tintQuiver);')
c_eval('positions? = xyz?+R?.resample(xyz?);')
c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?); hold(hca,''off'');')
hca.XLabel.String = 'N';
hca.YLabel.String = '-M';
hca.ZLabel.String = 'L';

hca = h2(4);
tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
c_eval('tts? = double((mvaVe?.tlim(tintQuiver).time.ttns-mvaVe?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
c_eval('R? = mvaR?-mvaR0;')
velocity = [0 0 30];
c_eval('xyz? = irf.ts_vec_xyz(mvaVe?.tlim(tintQuiver).time,tts?*velocity);')
c_eval('quivers? = mvaVe?.tlim(tintQuiver);')
c_eval('positions? = xyz?+R?.resample(xyz?);')
c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?,mms_colors(''?'')); hold(hca,''off'');')
hca.XLabel.String = 'N';
hca.YLabel.String = '-M';
hca.ZLabel.String = 'L';
end
%% Example 3 (Ohm's law, 4 sc): Plot 2x sc positions + 4x projection, ions 

tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:35.00Z');
%tt = irf.tint('2015-10-16T10:33:26.350Z',0.1);
%tt = irf.tint('2015-10-16T10:33:30.440Z',0.1)+(-0.18);
%tt = irf.tint('2015-10-16T10:33:27.920Z',0.1);
%tt=tt(1);

times = disDist1.tlim(tintZoom).time;


for it = 1:times.length;
  time = times(it);
  tint = time;
  tt = time;
  
  
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,tt.epochUnix','green');

  vectors = {'B',{mvaB1.resample(tt(1)),mvaB2.resample(tt(1)),mvaB3.resample(tt(1)),mvaB4.resample(tt(1))},2;...
             'Vi',{mvaVi1.resample(tt(1)),mvaVi2.resample(tt(1)),mvaVi3.resample(tt(1)),mvaVi4.resample(tt(1))},70};        

  isub = 1;         
  hca = h2(isub); isub = isub + 1; 
  lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 -1 0])
  axis(hca,'square')

  hca = h2(isub); isub = isub + 1; 
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 0 1])
  axis(hca,'square')

  %for ic = 1:4;
  %hca = h2(isub); isub = isub + 1; 
  %mms.plot_projection(hca,desDist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  vlim = 1*1e3;
    elevlim = 20;
    strCMap = 'jet';
    %energies =  [30 220];
    projclim = [3 9];
    palim = [1e-3 1e6];
    skymapEnergy = [65 278];


  for ic = 1:4;
      c_eval('dist = disDist?;',ic)

      c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
      hatVe0 = double(irf_norm(Ve0));    
      c_eval('Vi0 = gseVi?.resample(time).data;',ic); 
      hatVi0 = double(irf_norm(Vi0));    

      % Get mean magnetic field direction
      c_eval('B0 = gseB?.resample(time).data;',ic); 
      hatB0 = double(irf_norm(B0));
      c_eval('E0 = gseE?.resample(time).data;',ic); 
      hatE0 = double(irf_norm(E0));
      hatExB0 = cross(hatE0,hatB0);

      vectors = {hatB0,'B';hatE0,'E';hatVi0,'V_i';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

      % Projection coordinate system
      if 0
        x = hatB0;
        y = hatExB0;
        z = cross(x,y);
      elseif 0
        x = [1 0 0];
        y = [0 1 0];
        z = [0 0 1];        
      else        
        x = -N;
        y = L;
        z = -M;         
      end



      %isub = ic;

      if 0 % Plot psd 0 90 180
        hca = h2(isub); isub = isub + 1;
        psdtint = tint;%+[-1 1]*0.03;
        c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
        TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
        %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
        %hca.Title.String = irf_ssub('C?',ic);
        hca.Title.String = TitleStr;
      end

      if 0 % Plot skymap for a given energy
        hca = h2(isub); isub = isub + 1;      
        c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
        %hca.Title.String = hca.Title.String{2};

        % Plot skymap for a given energy
        hca = h2(isub); isub = isub + 1;      
        c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
        %hca.Title.String = hca.Title.String{2};
      end

      % Plot project ion onto a plane
      if 0
        hca = h2(isub); isub = isub + 1; 
        mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
        hca.Title.String = '';
        colormap(hca,strCMap)
      end

      hca = h2(isub); isub = isub + 1; 
       
      vlabel = {'-N','L','M'};
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',vlabel);
      titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
      hca.Title.String = titleStr;
      colormap(hca,strCMap)

      if 0
        hca = h2(isub); isub = isub + 1; 
        mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
        hca.Title.String = '';
        colormap(hca,strCMap)
      end
    end
    pause(1)
    cn.print([irf_ssub('ohm_sc_pos_4sc_ions',1) irf_time(time,'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end
  
if 0
hca = h2(3);
tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
c_eval('tts? = double((mvaB?.tlim(tintQuiver).time.ttns-mvaB?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
c_eval('R? = mvaR?-mvaR0;')
velocity = [0 0 30];
c_eval('xyz? = irf.ts_vec_xyz(mvaB?.tlim(tintQuiver).time,tts?*velocity);')
c_eval('quivers? = mvaB?.tlim(tintQuiver);')
c_eval('positions? = xyz?+R?.resample(xyz?);')
c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?); hold(hca,''off'');')
hca.XLabel.String = 'N';
hca.YLabel.String = '-M';
hca.ZLabel.String = 'L';

hca = h2(4);
tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
c_eval('tts? = double((mvaVe?.tlim(tintQuiver).time.ttns-mvaVe?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
c_eval('R? = mvaR?-mvaR0;')
velocity = [0 0 30];
c_eval('xyz? = irf.ts_vec_xyz(mvaVe?.tlim(tintQuiver).time,tts?*velocity);')
c_eval('quivers? = mvaVe?.tlim(tintQuiver);')
c_eval('positions? = xyz?+R?.resample(xyz?);')
c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?,mms_colors(''?'')); hold(hca,''off'');')
hca.XLabel.String = 'N';
hca.YLabel.String = '-M';
hca.ZLabel.String = 'L';
end
