function [Eper Epar] = cn_evinkb(diE,gseB)
%
% Returns diE that is always perpendicular to B or Bav
%
% E in ISR2
% bhat in GSE
% M transforms GSE to field aligned or vice versa
%

% get b in isr2-coordinates
diB=c_coord_trans('gse','isr2',gseB,'cl_id',3);

% get angle between B and spin plane
%cosalfa=sqrt(diB(:,2).^2+diB(:,3).^2);

%pardir=diB(:,1:2);
%parhat=pardir./sqrt(pardir(:,1).^2+pardir(:,2).^2);

% get isr2-coordinates of direction perpendicular to B and in the spin
% plane
if 0
for k=1:size(gseB);
    perdir=(cn_rotate(diB(k,2:3),90));
    perhat(k)=perdir/sqrt(perdir(1)^2+perdir(2)^2)
end
end

Eper=diE(:,1:2);
Epar=diE(:,1:2);
size(Eper);

for k=1:size(diE,1)
    parhat=diB(k,2:3)/sqrt(diB(k,2)^2+diB(k,3)^2);
    perdir=cn_rotate(diB(k,2:3),90);
    perhat=perdir/sqrt(perdir(1)^2+perdir(2)^2);
    cosalfa=sqrt(diB(k,2).^2+diB(k,3).^2)/sqrt(diB(k,2).^2+diB(k,3).^2+diB(k,4)^2);
    acosd(cosalfa);
    
    Eper(k,2)=(perhat*(diE(k,2:3))')';
    
    Epar(k,2)=parhat*(diE(k,2:3))'/cosalfa;
end

