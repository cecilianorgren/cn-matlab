function csc1 = tred_locate(xlim,ylim,zlim,Position)

[c1,sc1]=find(Position(:,:,1)>xlim(1));
[c2,sc2]=find(Position(:,:,1)<xlim(2));
csc1=intersect([c1,sc1],[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,2)>ylim(1));
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,2)<ylim(2));
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,3)>zlim(1));
csc1=intersect(csc1,[c2,sc2],'rows');
[c2,sc2]=find(Position(:,:,3)<zlim(2));
csc1=intersect(csc1,[c2,sc2],'rows');
