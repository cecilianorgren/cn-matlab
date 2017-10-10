function out = parseunderline(in)
% PARSEUNDERLINE Replaces '_' by '\_' for tex interpretation of strings.

indUnderline = strfind(in,'_');
cellOut = cellstr(tocolumn(in));
cellOut(indUnderline) = {'\_'};
matOut = char(cellOut');
spaceOut=reshape(matOut',1,numel(matOut));
indSpace = strfind(spaceOut,' ');
out = spaceOut;
out(indSpace) = [];
