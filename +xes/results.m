% saves results from xes1 simulations
% 'event' is the same as in xes.events
event = 'R0.20 S0.50';

switch event
    case 'R0.20 S0.50' % Daniels new MP event, slow ESWs, include in paper
        inputFile = 'event_case_r020s050.inp';        
        UEt = [2.475 8.046 22.9 36.1 48.6899 73.86];
        UE = [2.00534e-9 2.3767e-9 3.01325e-9 3.50395e-9 3.90181e-9 4.2068e-9];
    case 'R0.25 S0.50' % 
    case 'R0.30 S0.50' % 
        
    case 'R0.25 S0.40' % 
    case 'R0.25 S0.45' % 
        inputFile = 'event_case_r025s045.inp';        
        UEt = [15.26 64.70 92.57 132.04 186.16 297.24];
        UE = [2.72109e-9 5.20809e-9 5.22737e-9 7.25166e-9 8.77471e-9 1.17822e-8];
    case 'R0.20 S0.40' % 
        
    case 'R0.20 S0.35' % 
        
    otherwise
        disp('Did not recognize case name.')
        doInput = 0;
end

nt = numel(UEt);
for it = 1:(nt-1)
    growthrate(it) = log(UE(it+1)/UE(it))/(UEt(it+1)-UEt(it))/2;
    t(it) = (UEt(it+1)+UEt(it))/2;
end

subplot(2,1,1)
plot(t*sqrt(1/1836),growthrate*sqrt(1836),'-*')

subplot(2,1,2)
plot(UEt*sqrt(1/1836),UE,'-*')