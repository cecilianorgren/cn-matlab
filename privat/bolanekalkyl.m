loan = 1e6:1e5:3e6;
loan = 1e6;
interest = 0.04;
monthly_payment = 20000;

month = monthly_outcome(loan,interest,monthly_payment);
year = yearly_outcome(month);

% Run analysis for different loan amounts, monthly payments, and interests
loan = 1e6:1e5:3e6;
interest = 0.02:0.01:0.06;
monthly_payment = 15000:1000:20000;

for il = 1:numel(loan)
  for ii = 1:numel(interest) 
    for im = 1:numel(monthly_payment)
      month = monthly_outcome(loan(il),interest(ii),monthly_payment(im));
      outcome{il,ii,im} = yearly_outcome(month);
    end
  end
end

nrows = 3;
ncols = 3;
ipanel = 1;
for irow = 1:nrows
  for icol = 1:ncols
    ipanel = ipanel + 1;
    h(ipanel) = subplot(ncols,nrows,ipanel);
  end
end

isub = 1;




function out = monthly_outcome(loan,loan_interest,monthly_payment)
  
  months = [];
  loan_in = [];
  loan_out = [];
  interest = [];
  payoff = [];

  month = 0;
  while loan > 0
    month = month + 1;
    months(end+1) = month;
    loan_in(end+1) = loan;    
    
    
    monthly_interest = loan*loan_interest/12;
    monthly_loan_payoff = monthly_payment - monthly_interest;
    if monthly_loan_payoff > loan
      monthly_loan_payoff = loan;
    end
    loan = loan - monthly_loan_payoff;   

    loan_out(end+1) = loan;
    interest(end+1) = monthly_interest;
    payoff(end+1) = monthly_loan_payoff;
  end  

  out.month = months;
  out.loan_in = loan_in;
  out.loan_out = loan_out;
  out.interest = interest;
  out.payoff = payoff;
end

function out = yearly_outcome(in)
  fields = fieldnames(in);
  out.year = 1:ceil(numel(in.(fields{1}))/12);
  size_pad = 12-mod(numel(in.(fields{1})),12);

  tmp = [in.loan_in, zeros(1,size_pad)];
  out.loan_in = tmp(1:12:end);
  tmp = [in.loan_out, zeros(1,size_pad)];
  out.loan_out = tmp(12:12:end);
  
  for fields = ["interest","payoff"];
    tmp = [in.(fields), zeros(1,size_pad)];
    tmp = reshape(tmp,12,numel(tmp)/12);
    tmp = sum(tmp,1);
    out.(fields) = tmp;
  end
end