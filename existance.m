Ln=20*10/2.5

ro_p=427;
%%
c_eval('Neloc?=cn_toepoch(t1,t2,Ne?);',3:4);

dNe=abs(Neloc3(fix(end/2))-Neloc4(fix(end/2)));
Ln=0.1*10/dNe
%% 'figure
figure;irf_plot({Ne3,Ne4},'comp')
%%
PNeloc3
c_eval('PNeloc?=cn_toepoch(t1,t2,PNe?);',3:4);

dNe=abs(PNeloc3(fix(end/2))-PNeloc4(fix(end/2)));
Ln=0.08*10/dNe

%%
c_eval('PNeloc?=cn_toepoch(t1,t2,PNe?);',3:4);

dNe=abs(PNeloc3(end,2)-PNeloc4(end,2));
Ln=0.08*10/dNe
%%
Ln/ro_p
1863^(1/4)
1863^(1/2)
%%
figure;irf_plot({PNe3,PNe4},'comp')