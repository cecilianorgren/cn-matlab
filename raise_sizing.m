sb = 0.5;
bb = 1;
raise = 2.5*bb;
reraise = 2*raise;
call = reraise - raise;
pot = sb+bb+raise+reraise+call;

potodds = reraise/pot;

%%
sb = 0.5;
bb = 1;
raise = @(raisefactor) raisefactor*bb;
reraise = @(raisefactor,reraisefactor) reraisefactor.*raise(raisefactor);
call = @(raisefactor,reraisefactor) reraise(raisefactor,reraisefactor) - raise(raisefactor);
pot = @(raisefactor,reraisefactor) sb+bb+raise(raisefactor)+reraise(raisefactor,reraisefactor)+call(raisefactor,reraisefactor);

potodds = @(raisefactor,reraisefactor) reraise(raisefactor,reraisefactor)./pot(raisefactor,reraisefactor);

raisefactor_ = linspace(1,4,101);
reraisefactor_ = linspace(1,4,100);
[R,RER] = meshgrid(raisefactor_,reraisefactor_);

PODDS = potodds(R,RER);
hca = subplot(1,1,1);
[c,h] = contourf(hca,R,RER,PODDS);
clabel(c,h)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'My raise/(pot if opponent calls)';
shading(hca,'flat')
hca.XLabel.String = 'Opponent raise (bb)';
hca.YLabel.String = 'My reraise (op raise)';