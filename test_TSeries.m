% Test script for TSeries.mtimes

Matrix = irf.ts_tensor_xyz(gsePe1.time,gsePe1.data);
Vector = irf.ts_vec_xyz(gseB1.time,gseB1.data);
ScalarSingle = irf.ts_scalar(ne1.time,ne1.data);
ScalarVector = irf.ts_scalar(gseB1.time,gseB1.data);
ScalarMatrix = irf.ts_scalar(ePDist1.time,ePDist1.data);

L = [0.1357 0.1669 0.9766]; M = [0.3104 -0.9433 0.1180]; N = [0.9409 0.2871 -0.1798];

lmnVector1 = irf.ts_vec_xyz(Vector.time,[Vector.dot(L).data Vector.dot(M).data Vector.dot(N).data]);
lmnVector2 = Vector*[L' M' N'];

lmnMatrix1 = mms.rotate_tensor(Matrix,'rot',L,M,N); lmnMatrix1 = irf.ts_tensor_xyz(lmnMatrix1.time,lmnMatrix1.data);
lmnMatrix2 = [L; M; N]*Matrix*[L; M; N]';

%% Scalars
figure(3);
h = irf_plot(2);

isub = 1;

hca = h(isub); isub = isub + 1;
irf_plot(hca,{ScalarSingle*1.1,1.2*ScalarSingle,ScalarSingle*0,-ScalarSingle,ScalarSingle*-1.2},'comp')
hca.YLabel.String = {'Scalar, single'};

hca = h(isub); isub = isub + 1;
irf_plot(hca,{ScalarVector*[1 1 0]})
hca.YLabel.String = {'Scalar, vector'};

%% Vectors
figure(1);
h = irf_plot(2);

isub = 1;
hca = h(isub); isub = isub + 1;
irf_plot(hca,{lmnVector1.x,lmnVector1.y,lmnVector1.z},'comp')
hca.YLabel.String = {'Vector','LMN using dot'};

hca = h(isub); isub = isub + 1;
irf_plot(hca,{lmnVector2.x,lmnVector2.y,lmnVector2.z},'comp')
hca.YLabel.String = {'Vector','LMN using mtimes'};

if 0
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{lmnVector1.x-lmnVector2.x,lmnVector1.y-lmnVector2.y,lmnVector1.z-lmnVector2.z},'comp')
  hca.YLabel.String = {'Vector','difference'};
end

%% Matrices: Tensor order = 2 
figure(2);
h = irf_plot(4);

isub = 1;
hca = h(isub); isub = isub + 1;
irf_plot(hca,{lmnMatrix1.xx,lmnMatrix1.yy,lmnMatrix1.zz,lmnMatrix1.xy,lmnMatrix1.xz,lmnMatrix1.yz},'comp')
hca.YLabel.String = {'Tensor matrix','mms.rotate\_tensor'};

hca = h(isub); isub = isub + 1;
irf_plot(hca,{lmnMatrix2.xx,lmnMatrix2.yy,lmnMatrix2.zz,lmnMatrix2.xy,lmnMatrix2.xz,lmnMatrix2.yz},'comp')
hca.YLabel.String = {'Tensor matrix','R*T*R'''};

if 1 % Diagonal terms
  hca = h(isub); isub = isub + 1;  
  lines1 = irf_plot(hca,{lmnMatrix1.xx,lmnMatrix1.yy,lmnMatrix1.xz},'comp');
  hold(hca,'on')
  lines2 = irf_plot(hca,{lmnMatrix2.xx,lmnMatrix2.yy,lmnMatrix2.xz},'comp');  
  hold(hca,'on')
  
  for ii = 1:3;
    lines2.Children(ii+3).Color = lines2.Children(ii).Color;
    lines2.Children(ii+3).LineStyle = '--';
  end
  hca.YLabel.String = {'Diagonal terms'};
  irf_legend(hca,{'mms.rotate\_tensor','R*T*R'''},[0.02 0.98])
end

if 1 % Off-diagonal terms
  hca = h(isub); isub = isub + 1;
  lines1 = irf_plot(hca,{lmnMatrix1.xy,lmnMatrix1.xz,lmnMatrix1.yz},'comp');
  hold(hca,'on')
  lines2 = irf_plot(hca,{lmnMatrix2.xy,lmnMatrix2.xz,lmnMatrix2.yz},'comp');  
  hold(hca,'on')
  
  for ii = 1:3;
    lines2.Children(ii+3).Color = lines2.Children(ii).Color;
    lines2.Children(ii+3).LineStyle = '--';
  end
  hca.YLabel.String = {'Off-diagonal terms'};
  irf_legend(hca,{'mms.rotate\_tensor','R*T*R'''},[0.02 0.98])
end

%%
figure(3)

plot(lmnMatrix1.xx.data,lmnMatrix2.xx.data,'.',...
     lmnMatrix1.yy.data,lmnMatrix2.yy.data,'.',...
     lmnMatrix1.zz.data,lmnMatrix2.zz.data,'.',...
     lmnMatrix1.xy.data,lmnMatrix2.xy.data,'.',...
     lmnMatrix1.xz.data,lmnMatrix2.xz.data,'.',...
     lmnMatrix1.yz.data,lmnMatrix2.yz.data,'.',...
     [0 0.2],[0 0.2],'k--')

%% Partial tensors: implementing plus(), minus(), mtimes(), and rmdivide()




