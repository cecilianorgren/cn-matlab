
function out_f=Zplasma_d(in_f)
%out_f=i*sqrt(pi)*exp(-in_f.*in_f).*(-2*in_f).*erfc_complex(-i*in_f)-2;

% if(abs(in_f)<30)
% 
% out_f=-2*(1+(in_f.*Zplasma(in_f)));
% 
% else
%    
% %   out_f= - 2*i*in_f*sqrt(pi)*exp(-(in_f*in_f))
%     out_f= (1+(1.5/(in_f*in_f)))/(in_f*in_f);
%    
% end





%  if(abs(in_f)<20)
%  
%  out_f=-2*(1+(in_f.*Zplasma(in_f)));
%  
%  else
%     
% %%%   out_f= - 2*i*in_f*sqrt(pi)*exp(-(in_f*in_f))
%      out_f= (1+(1.5/(in_f*in_f)))/(in_f*in_f);
%     
%  end

out_f=-2*(1+(in_f.*Zplasma(in_f)));


end
