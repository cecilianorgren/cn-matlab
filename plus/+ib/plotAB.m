% Will compare the data that has been had from ib.A and ib.B respectively.
% the following .mat files has been created
% mA.mat            Atwo4
% mEDSI.mat         
% mEDSIf.mat        
% mEFW.mat          
% mEFWR.mat         
% mEFWburst.mat     
% mEFWburstR.mat    
% mER.mat           
% mP.mat            
% mPR.mat           
% mEFWburstTM.mat (only B (caa_get_burst))

direct = {'/Users/Cecilia/Research/Faz/data/A',...
          '/Users/Cecilia/Research/Faz/data/B'};
         % '/Users/Cecilia/Research/Faz/data/C'};
way  = {'A','B'};
mats = {'mA.mat',...
        'mEDSI.mat',...
        'mEDSIf.mat',...
        'mEFW.mat',...    
        'mEFWR.mat',...      
        %'mEFWburst.mat',...      
        'mEFWburstR.mat',...      
        'mER.mat',...      
        'mP.mat',...      
        'mPR.mat',...       
        'mEFWburstTM.mat'};
nMats = numel(mats);

if 1
fid = fopen('/Users/Cecilia/Research/Faz/data/AB/products_A.txt','w');

for ii = 1:nMats-1 % runs over different .mat
    clearvars -except ii mats nMats direct fid
    for jj = 1%:2 % runs over A and B method
        cd(direct{jj});
        load(mats{ii});    
        vars = whos; 
        nVars = numel(vars);        
        for kk = 1:nVars
            if max(vars(kk).size) > 13 % probably a time series
                % Add variable to list in order to later plot and compare 
                % them. Load with c_load or just use irf_plot directly.
                strVar = vars(kk).name;
                strMat = mats{ii};
                %fprintf(fid,[strVar strMat],'%s %*f %s \n',[strVar strMat]);  
                fprintf(fid,[sprintf('%-*s',30,strVar) strMat '\n']);  
            end
        end
    end
end
type '/Users/Cecilia/Research/Faz/data/AB/products.txt';
fclose(fid);
end

if 0    
    for ii = 1:nMats-1
        for jj = 1:2
            cd(direct{jj});
            load(mats{ii}); 
        end
    end
end
    

%% Make plots to compare
fid = fopen('/Users/Cecilia/Research/Faz/data/AB/products_A.txt','r');
list = textscan(fid,'%s%s');
nProd = numel(list{1});

for pp = 1:nProd
    var = char(list{1}(pp));
    matFile = list{2}{pp};
    h = irf_plot(3);
    
    % method A
    cd(direct{1}); 
    load(matFile); 
    eval(['Atemp = ' var ';']);
    clear(var);
    irf_plot(h(1),Atemp);   
    ylabel(h(1),'A')
    
    % method B
    cd(direct{2}); 
    load(matFile); 
    eval(['Btemp = ' var ';']);
    clear(var);
    irf_plot(h(2),Btemp);
    ylabel(h(2),'B   (caa\_get\_bursts)')
    
    % difference between A and B
    ABtemp = irf_add(1,Atemp,-1,Btemp);
    irf_plot(h(3),ABtemp);
    ylabel(h(3),'A-B')
    
    % title and print
    title(h(1),[var '  from  ' matFile])
    cn.print(var);
end





