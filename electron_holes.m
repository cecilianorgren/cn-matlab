x=b_hat;
y=n_hat;
z=irf_cross(b_)
%%
% Reconstructing Etot assuming it is parallel field
c_eval('b?=cn_toepoch(t1,t2,gsmB?);',3:4)
c_eval('[gsmEt? eb_ang?]=irf_edb(gsmE?(:,1:3),b?,90,''Epar'');',3:4);
c_eval('gsmBw?=irf_filt(gsmB?,1,0,450,5);',3:4);
%c_eval('gsmEw?=irf_filt(gsmEt?,1,0,450,5);',3:4);

% Transform to field-aligned coordinates
z=[b3(:,1) b3(:,2:4)./repmat(b3(:,5),size(b3,1),size(b3,2))];

% Plot E and B
h=irf_plot(2);
irf_plot(h(1),gsmEt3)
irf_plot(h(2),gsmBw3)
%%
c_eval('b_hat?=[diB?(:,1) diB?(:,2:4)./repmat(diB?(:,5),1,3)];',3:4)
c_eval('Epar?=irf_dot(Etot?,b_hat?);',3:4)

c_eval('bxn?=irf_cross(b_hat?,[b_hat?(:,1) repmat([0 0 1],size(b_hat?,1),1)]);',3:4)
c_eval('Eper?=irf_dot(Etot?,bxn?);',3:4)