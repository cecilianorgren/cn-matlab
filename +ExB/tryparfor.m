% try parfor
nl = 10;
rn = zeros(nl,100);
parfor ii = 1:nl
    rn(ii,:) = rand(100,1);
end
    
%%
nl = 100;
nk = 10;
rn = zeros(nl,nk);
num = ones(nl);
%rn2 = rn*0;
si = 1:fix(nl/nk):101;
vect = zeros(nl,nk);
fillin = zeros(100,1);
parfor kk = 1:nk;
    parfillin = zeros(100,1);
    nfor = si(kk):(si(kk+1)-1);
    for ll = nfor
        new = round(100*rand(1,1));        
        %disp(num2str([ll kk]))
        parfillin(new) = parfillin(new) + 1;        
    %rn(ll,kk) = ll + num();
    end
    fillin = fillin + parfillin;
end