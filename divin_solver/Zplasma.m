
function out_f=Zplasma(in_f)

%out_f=(i*sqrt(pi)*exp(-in_f.*in_f).*erfc_complex(-i*in_f)).*(imag(in_f)>=0)+(conj(i*sqrt(pi)*exp(-conj(in_f).*conj(in_f)).*erfc_complex(-i*conj(in_f)))+2*i*sqrt(pi)*exp(-conj(in_f).*conj(in_f) ) ).*(imag(in_f)<0);
%out_f=i*sqrt(pi)*exp(-in_f.*in_f).*erfc_complex(-i*in_f);

%  if(abs(in_f)<20)
%    out_f = i*sqrt(pi).*exp(-in_f.*in_f).*mfun('erfc',-i*in_f);
% 
%  else
%      out_f = (-( 1+(0.5./(in_f.^2))+(0.75./(in_f.^4)))./in_f)+(1i*sqrt(pi)*exp(-in_f.^2)  );
%  end




   out_f =  i*sqrt(pi)*Faddeeva(in_f,64);
  %out_f =  i*sqrt(pi)*Faddeeva_another_version(in_f);


%if (isinf(exp(-in_f.*in_f)) || isinf(mfun('erfc',-i*in_f))  ) out_f=0; end
%if (isnan(out_f)) out_f=0;
    
end
