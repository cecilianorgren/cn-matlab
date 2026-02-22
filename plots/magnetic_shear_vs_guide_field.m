B0 = 10;
BG = (0:10)';
nBG = numel(BG);

B1 = [repmat(B0,nBG,1) BG ];
B2 = [repmat(-B0,nBG,1) BG];

B1abs = sqrt(sum(B1.^2,2));
B2abs = sqrt(sum(B2.^2,2));

B1 = B1./repmat(B1abs,1,2);
B2 = B2./repmat(B2abs,1,2);

dotb1b2 = sum(B1.*B2,2);
shear = acosd(dotb1b2);

plot(BG/B0,shear)