function [wr_out,wi_out,wr_allout,wi_allout,wrin_out,wiin_out,wrin_allout,wiin_allout,param_out] = beam_solver(Te1,Te2,Ti,R,S,k,n,doPlot,TolFun)
% Disperison function solver.
%   [wr,wi,parameters] = beam_solver(Te1,Te2,Ti,R,S,k);
%
%
%
%
%

TolFun
if isempty(doPlot)
    doPlot = 0;
end

% input temperatures is normalized to Te1, so turning it back now.
Te2=Te1*Te2;
Ti=Te1*Ti;

nR = numel(R);
nv = numel(S);
nk = numel(k);
nTe1 = numel(Te1);
nTe2 = numel(Te2);
nTi = numel(Ti);

disp('Running for:')
disp(['k = ' num2str(k)])
disp(['R = ' num2str(R)])
disp(['S = ' num2str(S)])
disp(['Te1 = ' num2str(Te1)])
disp(['Te2 = ' num2str(Te2)])
disp(['Ti = ' num2str(Ti)])

B = 25;
no = 0;

%units = irf_units;
e = 1.6022e-19; %units.e;
eps0 = 8.8542e-12; %units.eps0;
mp = 1.6726e-27; %units.mp;
me = 9.1094e-31; %units.me;
c = 299792458; %units.c;

omega_pi = sqrt(n*1e6*e^2/eps0/mp); % rad/s
omega_pe = sqrt(n*1e6*e^2/eps0/me); % rad/s
omega_pe1 = sqrt(n*1e6*e^2/eps0/me)*sqrt(1-R); % rad/s
omega_pe2 = sqrt(n*1e6*e^2/eps0/me)*sqrt(R); % rad/s

vthe1 = c*sqrt(1-1./(Te1.*e./(me*c^2)+1).^2);
vthe2 = c*sqrt(1-1./(Te2.*e./(me*c^2)+1).^2);
vthi = c*sqrt(1-1./(Ti.*e./(mp*c^2)+1).^2);

vde1 = 0; % m/s
vde2 = S*vthe1; % m/s
vdi = 0; % m/s

lamD1 = (vthe1')*(1./omega_pe1)/sqrt(2); % m
lamD2 = (vthe2')*(1./omega_pe2)/sqrt(2); % m

lamD = (vthe1')*(1./omega_pe)/sqrt(2); % m
%l_norm = v_norm/w_norm/sqrt(2);
% omega_bune = omega_pe^(1/3)*omega_pi^(2/3)*16^(-1/3);

normal_frequency = omega_pe; % total electron plasma frequency
normal_length = lamD; % Vte./Wpe/sqrt(2); % total density, background temperature
normal_frequency_string = '\omega_{pe}';
normal_length_string = '\lambda_D';

param_out.B = B;
param_out.n = n;
param_out.ope = omega_pe;
param_out.opi = omega_pi;
param_out.l_Debye = lamD;
param_out.l_Debye1 = lamD1;
param_out.l_Debye2 = lamD2;
param_out.vte1 = vthe1;
param_out.vte2 = vthe2;
param_out.vti = vthi;
param_out.wnorm = normal_frequency;
param_out.lnorm = normal_length;

if 0
wr_out=0;
wi_out=0;
return
end
if doPlot % set up figure for drawing results 
    
    figure(10)    
    for nPlot = 1:2     
        h(nPlot) = subplot(1,2,nPlot); hold(h(nPlot),'on');        
        xlabel(h(nPlot),['k' normal_length_string],'Fontsize',14);
        ylabel(h(nPlot),['\omega/'  normal_frequency_string],'Fontsize',14);    
        set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
    end
    title(h(1),'\Re (\omega)')
    title(h(2),'\Im (\omega)')
end

nTryLow = 1;
nTryHigh = 0;

rights = zeros(1,nk);

% initialize matrices
%x_real_store = zeros(nk,nv,nR,nTe2,nTi);
%x_imag_store = zeros(nk,nv,nR,nTe2,nTi);

% Starting guesses, can find it from solver_BunemanRandom
x_real_store(1,:,:,:,:) = ones(1,nv,nR,nTe2,nTi,1)*0.0001*normal_frequency;
x_imag_store(1,:,:,:,:) = ones(1,nv,nR,nTe2,nTi,1)*0.0005*normal_frequency;
x_real_storeall(1,:,:,:,:) = nan;
x_imag_storeall(1,:,:,:,:) = nan;

x_real_in_storeall = nan;
x_imag_in_storeall = nan;
%
for iTi = 1:nTi;
    for iTe2 = 1:nTe2;              
        for iR = 1:nR % Different beam densities            
            for iv=1:nv;          
                disp(['Ti=' num2str(Ti(iTi)) ' Te2=' num2str(Te2(iTe2)) ' R=' num2str(R(iR)) ' S=' num2str(S(iv))])
                surfOK = 0; % has not yet found a positive growth rate surface
                wimax = 1;
                for ik=1:nk % For each k*lD    
                    fv = [omega_pi,omega_pe*sqrt(1-R(iR)),omega_pe*sqrt(R(iR)),vthi(iTi),vthe1,vthe2(iTe2),vdi,vde1,vthe1*S(iv),lamD];
                    af = @(temp) beam_solver_disprel(temp,k(ik),fv);
                    if k(ik)<0.13; TolFun=1e-8; elseif k<0.5; TolFun=1e-13; else TolFun = 1e-14; end
                    
                    initialGuess = 5;
                    switch initialGuess
                        case 1 % old guessing structure
                            % make appropriate initial guess
                            if ik == 1; 
                                win = [x_real_store(1,iv,iR,iTe2,iTi) x_imag_store(1,iv,iR,iTe2,iTi)];
                            else
                                if S(iv) > 1 %&& x_real_store(ik,iv-1,iR,iTe2,iTi)<0
                                    win = [x_real_store(ik,iv-1,iR,iTe2,iTi) x_imag_store(ik,iv-1,iR,iTe2,iTi)];
                                else
                                    if x(1) < 0 || x(1) > 0.95*normal_frequency;
                                        win = [0.001 0.0005]*normal_frequency;
                                    else
                                        if iv == 1

                                            %iv
                                            win = [x_real_store(ik-1,1,iR,iTe2,iTi) x_imag_store(ik-1,1,iR,iTe2,iTi)];
                                        else
                                            %ik 
                                            %iv                                    
                                            win = [x_real_store(ik-1,iv-1,iR,iTe2,iTi) x_imag_store(ik-1,iv-1,iR,iTe2,iTi)];;
                                        end
                                    end
                                end
                            end


                            %win = [x_real_store(kend,iv,iR,iTe2,iTi) x_imag_store(kend,iv,iR,iTe2,iTi)];


                            [x,FVAL,EXITFLAG] = fsolve(af,win,optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));
                            %if(x(1)>0 && x(1)<normal_frequency*0.8)                        
                                x_real_store(ik,iv,iR,iTe2,iTi) = x(1);
                                x_imag_store(ik,iv,iR,iTe2,iTi) = x(2);
                                %disp(['xin = ' num2str(win) ', xout = ' num2str(x)]);% ', aftry = ' num2str(aftry)]) 
                            %end                                        
                        case 2 % new guessing structure
                            % make appropriate initial guess
                            if ik == 1; 
                                win = [x_real_store(1,iv,iR,iTe2,iTi) x_imag_store(1,iv,iR,iTe2,iTi)];
                            else
                                if isOk % found positive wi last round, use this as initial guess for next k
                                    win = [x_real_store(ik-1,iv,iR,iTe2,iTi) x_imag_store(ik-1,iv,iR,iTe2,iTi)];
                                else % did not find wi>0
                                    % use last S-solution for same k as initial guess
                                    if iv>1
                                        win = [x_real_store(ik,iv-1,iR,iTe2,iTi) x_imag_store(ik,iv-1,iR,iTe2,iTi)];
                                    else % use the very first guess
                                        win = [x_real_store(1,iv,iR,iTe2,iTi) x_imag_store(1,iv,iR,iTe2,iTi)];
                                    end
                                end                        
                            end

                            [x,FVAL,EXITFLAG] = fsolve(af,win,optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));                                      
                            x_real_store(ik,iv,iR,iTe2,iTi) = x(1);
                            x_imag_store(ik,iv,iR,iTe2,iTi) = x(2);   

                            if ik==1
                                isOk = 1;
                            else % ik ~= 1 
                                if x(2)>0; isOk = 1; 
                                elseif any(x_imag_store(1:ik,iv,iR,iTe2,iTi))>0 && ik>3
                                    % if any point on the surface was positive,
                                    % continue on it
                                    isOk = 1;
                                elseif abs(x(2)-x_imag_store(ik-1,iv,iR,iTe2,iTi))/abs(x_imag_store(ik-1,iv,iR,iTe2,iTi))<0.5
                                    isOk = 1;
                                else
                                    isOk = 0;
                                end                    
                            end     
                            %if S(iv)< 0.5 isOk = 1; end                            
                        case 3 % newer
                            if ik == 1; 
                                win = [x_real_store(1,iv,iR,iTe2,iTi) x_imag_store(1,iv,iR,iTe2,iTi)];
                            else
                                if isOk % found positive wi last round, use this as initial guess for next k
                                    win = [x_real_store(ik-1,iv,iR,iTe2,iTi) x_imag_store(ik-1,iv,iR,iTe2,iTi)];
                                else % did not find wi>0
                                    % use last S-solution for same k as initial guess
                                    if iv>1
                                        win = [x_real_store(ik,iv-1,iR,iTe2,iTi) x_imag_store(ik,iv-1,iR,iTe2,iTi)];
                                    else % use the very first guess
                                        win = [x_real_store(1,iv,iR,iTe2,iTi) x_imag_store(1,iv,iR,iTe2,iTi)];
                                    end
                                end                        
                            end

                            [x,FVAL,EXITFLAG] = fsolve(af,win,optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));                                      
                            x_real_store(ik,iv,iR,iTe2,iTi) = x(1);
                            x_imag_store(ik,iv,iR,iTe2,iTi) = x(2);   

                            if ik==1
                                isOk = 1;
                            else % ik ~= 1 
                                if x(2)>0; isOk = 1; 
                                elseif any(x_imag_store(1:ik,iv,iR,iTe2,iTi))>0
                                    % if any point on the surface was positive,
                                    % continue on it
                                    isOk = 1;
                                elseif abs(x(2)-x_imag_store(ik-1,iv,iR,iTe2,iTi))/abs(x_imag_store(ik-1,iv,iR,iTe2,iTi))<0.5
                                    isOk = 1;
                                else
                                    isOk = 0;
                                end                    
                            end     
                            %if S(iv)< 0.5 isOk = 1; end
                        case 4 % newer
                            % make different guesses and store the one with
                            % highest imaginary frequency + the any surface 
                            % which has had a positive growth rate at any
                            % point
                            
                            clear winn xxr xxi xx 
                            % always have these guesses                                                     
                            winn{1}=[0.0001 0.0005]*normal_frequency; 
                            winn{2}=[0.01 0.005]*normal_frequency;  
                            if surfOK; 
                                winn{end+1} = [x_real_storeall(ik-1,iv,iR,iTe2,iTi,indexOK),...
                                               x_imag_storeall(ik-1,iv,iR,iTe2,iTi,indexOK)];
                                winn{end+1} = [x_real_storeall(ik-1,iv,iR,iTe2,iTi,indexOK),...
                                               x_imag_storeall(ik-1,iv,iR,iTe2,iTi,indexOK)*3];                                
                                winn{end+1} = [-0.5*normal_frequency,...
                                               x_imag_storeall(ik-1,iv,iR,iTe2,iTi,indexOK)];           
                            end
                            if ik == 1 % first k                                                                                                
                                if iv ~= 1 % not first S
                                    % same k, last S
                                    % highest growthrate
                                    winn{end+1} = [x_real_store(ik,iv-1,iR,iTe2,iTi,1) x_imag_store(ik,iv-1,iR,iTe2,iTi,1)];
                                end
                            else % not first k
                                % last point
                                winn{end+1} = [x_real_store(ik-1,1,iR,iTe2,iTi) x_imag_store(ik-1,1,iR,iTe2,iTi)];
                                if k > 2  
                                    if surfOK
                                        % interpolated from last point        
                                        dwr = x_real_storeall(ik-2,1,iR,iTe2,iTi,indexOK) - x_real_storeall(ik-1,1,iR,iTe2,iTi,indexOK);
                                        dwi = x_imag_storeall(ik-2,1,iR,iTe2,iTi,indexOK) - x_imag_storeall(ik-1,1,iR,iTe2,iTi,indexOK);
                                        winn{end+1} = [x_real_storeall(ik-1,1,iR,iTe2,iTi,indexOK)+dwr x_imag_storeall(ik-1,1,iR,iTe2,iTi,indexOK)+dwi];
                                    end
                                end
                                if iv == 1 % first S
                                    
                                elseif iv == 2 % second S
                                    winn{end+1} = [x_real_store(ik,iv-1,iR,iTe2,iTi) x_imag_store(ik,iv-1,iR,iTe2,iTi)];
                                else 
                                    if surfOK
                                        winn{end+1} = [x_real_storeall(ik,iv-1,iR,iTe2,iTi,indexOK),...
                                                       x_imag_storeall(ik,iv-1,iR,iTe2,iTi,indexOK)];
                                    end
                                    winn{end+1} = [x_real_store(ik,iv-1,iR,iTe2,iTi) x_imag_store(ik,iv-1,iR,iTe2,iTi)];
                                    winn{end+1} = [x_real_store(ik,iv-2,iR,iTe2,iTi) x_imag_store(ik,iv-2,iR,iTe2,iTi)];
                                end
                            end
                            for uu = 1:numel(winn)
                                x_real_in_storeall(ik,iv,iR,iTe2,iTi,uu) = winn{uu}(1);
                                x_imag_in_storeall(ik,iv,iR,iTe2,iTi,uu) = winn{uu}(2);
                            end
                            
                            if 1 % new selection, dont take in the all highest values
                                for ii = 1:numel(winn)  
                                    [x,FVAL,EXITFLAG] = fsolve(af,winn{ii},optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off'));
                                    %[x,FVAL,EXITFLAG] = lsqnonlin(af,winn{ii},0,10000,optimoptions('lsqnonlin','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off'));
                                    x_real_storeall(ik,iv,iR,iTe2,iTi,ii) = x(1);
                                    x_imag_storeall(ik,iv,iR,iTe2,iTi,ii) = x(2);
                                    
                                    % check if solver has found a positive
                                    % growth rate surface, if yes, continue
                                    % searching along this surface
                                    if x(2)>wimax; wimax = x(2); surfOK = 1; indexOK = ii; end
                                        
                               
                                    if x(2) < 10*normal_frequency % abnormally high, dont save
                                        xxr(ii) = x(1);
                                        xxi(ii) = x(2);
                                        xx{ii} = x;
                                    else % don't save if it's abnormally high growth rate
                                        xxr(ii) = NaN;
                                        xxi(ii) = NaN;
                                        xx{ii} = [xxr(ii) xxi(ii)];
                                    end
                                end                 
                                %try
                                %[isorted,isortind] = sort(xxi);
                                %rsorted = xxr(isortind);
                                
                                % expecially save solution with highest
                                % growth rate
                                maxind = find(xxi == max(xxi));
                                maxind = maxind(1);
                                x_real_store(ik,iv,iR,iTe2,iTi) = xxr(maxind);
                                x_imag_store(ik,iv,iR,iTe2,iTi) =  xxi(maxind);
                                
                                
                                %x_real_store(ik,iv,iR,iTe2,iTi,1) = rsorted(end);
                                %x_imag_store(ik,iv,iR,iTe2,iTi,1) = isorted(end);
                                %x_real_store(ik,iv,iR,iTe2,iTi,2) = rsorted(end-1);
                                %x_imag_store(ik,iv,iR,iTe2,iTi,2) = isorted(end-1);
                                win = winn{maxind};
                                %catch
                                %    1
                                %end
                            else % old selection
                                                                
                                                                
                                % additional guesses
                                for ii = 2:numel(winn)  
                                    [x,FVAL,EXITFLAG] = fsolve(af,winn{ii},optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));
                                    %if x(2)>xrtemp && x(2)
                                    xxr(ii) = x(1);
                                    xxi(ii) = x(2);
                                    xx{ii} = x;
                                end                 
                                try
                                maxind = find(xxi == max(xxi));

                                %maxind = maxind(1);

                                x_real_store(ik,iv,iR,iTe2,iTi) = xxr(maxind);
                                x_imag_store(ik,iv,iR,iTe2,iTi) = xxi(maxind);
                                win = winn{maxind};
                                catch
                                    1;
                                end
                                
                                
                            end
                        case 5 % take highest, but not abnormally high
                            % make different guesses and store the one with
                            % highest imaginary frequency + the any surface 
                            % which has had a positive growth rate at any
                            % point
                            
                            clear winn xxr xxi xx 
                            % always have these guesses                                                     
                            winn{1}=[0.0001 0.0005]*normal_frequency; 
                            winn{2}=[0.01 0.005]*normal_frequency;  
                            if 0;surfOK; 
                                winn{end+1} = [x_real_storeall(ik-1,iv,iR,iTe2,iTi,indexOK),...
                                               x_imag_storeall(ik-1,iv,iR,iTe2,iTi,indexOK)];        
                            end
                            if ik == 1 % first k                                                                                                
                                if iv ~= 1 % not first S
                                    % same k, last S
                                    % highest growthrate
                                    winn{end+1} = [x_real_store(ik,iv-1,iR,iTe2,iTi,1) x_imag_store(ik,iv-1,iR,iTe2,iTi,1)];
                                end
                            else % not first k
                                % last point
                                winn{end+1} = [x_real_store(ik-1,1,iR,iTe2,iTi) x_imag_store(ik-1,1,iR,iTe2,iTi)];                                                              
                            end
                            for uu = 1:numel(winn)
                                x_real_in_storeall(ik,iv,iR,iTe2,iTi,uu) = winn{uu}(1);
                                x_imag_in_storeall(ik,iv,iR,iTe2,iTi,uu) = winn{uu}(2);
                            end
                            
                            if 1 % new selection, dont take in the all highest values
                                for ii = 1:numel(winn)  
                                    [x,FVAL,EXITFLAG] = fsolve(af,winn{ii},optimoptions('fsolve','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off'));
                                    %[x,FVAL,EXITFLAG] = lsqnonlin(af,winn{ii},0,10000,optimoptions('lsqnonlin','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off'));
                                    x_real_storeall(ik,iv,iR,iTe2,iTi,ii) = x(1);
                                    x_imag_storeall(ik,iv,iR,iTe2,iTi,ii) = x(2);
                                    
                                    % check if solver has found a positive
                                    % growth rate surface, if yes, continue
                                    % searching along this surface
                                    if x(2)>wimax; wimax = x(2); surfOK = 1; indexOK = ii; end                                                                       
                                    if x(2) < 10*normal_frequency % abnormally high, dont save
                                        xxr(ii) = x(1);
                                        xxi(ii) = x(2);
                                        xx{ii} = x;
                                    else % don't save if it's abnormally high growth rate
                                        xxr(ii) = NaN;
                                        xxi(ii) = NaN;
                                        xx{ii} = [xxr(ii) xxi(ii)];
                                    end
                                end                 
                                
                                % expecially save solution with highest
                                % growth rate
                                maxind = find(xxi == max(xxi));
                                maxind = maxind(1);
                                x_real_store(ik,iv,iR,iTe2,iTi) = xxr(maxind);
                                x_imag_store(ik,iv,iR,iTe2,iTi) =  xxi(maxind);
                                
                                
                                win = winn{maxind};
                            end
                    end                                
                                                     
                win_real_store(ik,iv,iR,iTe2,iTi) = win(1); 
                win_imag_store(ik,iv,iR,iTe2,iTi) = win(2);                         
                if doPlot
                        if(EXITFLAG==1) % Found something, now plot            
                            plot(h(1),k(ik),x(1)/normal_frequency,'g+','Markersize',3);             
                            plot(h(2),k(ik),x(2)/normal_frequency,'g+','Markersize',3); hold on;
                            % only plot for reasonably omega_i > -0.6            

                            if(x(2)/normal_frequency>-0.6)
                                if ( abs(x(1)/normal_frequency)<0.8)                    
                                    plot(h(1),k(ik),x(1)/normal_frequency,'r+','Markersize',6); hold on; 
                                    plot(h(2),k(ik),x(2)/normal_frequency,'r+','Markersize',6); hold on;                
                                elseif 0;( abs(x(1)/normal_frequency)>0.8) % before it was 0.1                    
                                    plot(h(1),k(ik),x(1)/normal_frequency,'bd','Markersize',6); hold on;                    
                                    plot(h(2),k(ik),x(2)/normal_frequency,'bd','Markersize',6); hold on;    
                                end
                                if (x(2)/normal_frequency>0.0001)                    
                                    rights(ik)=rights(ik)+1;
                                    plot(h(1),k(ik),x(1)/normal_frequency,'ch','Markersize',6); hold on;                    
                                    plot(h(2),k(ik),x(2)/normal_frequency,'ch','Markersize',6); hold on; 
                                    %disp(['x_real=' num2str(x_real_store(end,iv,iR),'%.4f') ', x_imag=' num2str(x_imag_store(end,iv,iR),'%.4f')])
                                end
                            end                   
                            drawnow  
                        end                
                end
                size(x_imag_storeall);    
                end % end nk                        
            end % end nv            
        end % end nR
    end % end nTe2
end % end nTi

if 0
    set(h(1),'ylim',[0 1],'xlim',[0 2]);
    set(h(2),'ylim',[0 0.5],'xlim',[0 2])
end
if 0 % this can be done outside
    % make and store distributions
    ves = linspace(-10000,vthe1*3*1e-3,200);
    vis = linspace(-10000,vthe1*0.3*1e-3,200);
    fe_store = zeros(numel(ves),nv,nR);
    %fi_store = zeros(numel(vis),1,1);
    fi_store = cn.maxwellian(vis,Ti,n,0,'p',1); 

    for iR = 1:nR % Different beam densities    
        for iv=1:nv % For each beam drift speed        
            fe_store(:,iv,iR) = cn.maxwellian(ves,Te1,n*(1-R(iR)),0,'e',1)+cn.maxwellian(ves,60,n*R(iR),S(iv)*vthe1*1e-3,'e',1);

        end
    end
    %plot(ves,squeeze(fe_store(:,2,:)),vis,squeeze(fi_store(:,2,:))/10) 
end

wr_out = x_real_store/normal_frequency;
wi_out = x_imag_store/normal_frequency;
wr_allout = x_real_storeall/normal_frequency;
wi_allout = x_imag_storeall/normal_frequency;
wrin_out = win_real_store/normal_frequency;
wiin_out = win_imag_store/normal_frequency;
wrin_allout = x_real_in_storeall/normal_frequency;
wiin_allout = x_imag_in_storeall/normal_frequency;
