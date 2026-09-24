function varargout = gmm_bic_find_knee(Ks,BIC,threshold)

dBIC = -diff(BIC,1,2);
relativeDrop = dBIC(:,2:end) ./ dBIC(:,1:end-1);
relativeDropThreshold = 0.3;


nt = size(BIC,1);

BICreldr_K = zeros(nt,1);
BICreldr_BIC = zeros(nt,1);

BICredr_K_last = zeros(nt,1);
BICredr_BIC_last = zeros(nt,1);

for it = 1:nt

  % Option 4
  [idx] = find(relativeDrop(it,:)<threshold,1,'first');
  if isempty(idx), idx = 1; end
  BICredr_K(it) = Ks(idx+1);
  BICredr_BIC(it) = BIC(it,idx+1);

  % Option 5, also relative drop, but find the last one
  relativeDrop_pos = relativeDrop(it,:);
  relativeDrop_pos(relativeDrop_pos<0) = NaN;
  [idx] = find(relativeDrop_pos<threshold,1,'last');
  if isempty(idx)
    BICreldr_K(it) = NaN;
    BICreldr_BIC(it) = NaN;
  else
    BICreldr_K(it) = Ks(idx+1);
    BICreldr_BIC(it) = BIC(it,idx+1);
  end
end

