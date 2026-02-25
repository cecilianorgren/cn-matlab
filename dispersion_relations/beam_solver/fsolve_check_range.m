f = @(x) (x+0.5).^2+1;
fsolve(f,1+0.5i) % finds -0.5+1i
fsolve(f,1-1.5i) % finds -0.5-1i
%%
lsqnonlin(f,0.1+0.5i,-1,0,optimoptions('lsqnonlin','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off')) % should find -0.5+1i
sqrt(real(1).^2+imag(2).^2)
lsqnonlin(f,-0.7-1.5i,-1,0,optimoptions('lsqnonlin','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off')) % should find -0.5-1i
sqrt(real(1).^2+imag(2).^2)