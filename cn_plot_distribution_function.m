function cn_plot_distribution_function(varargin)
% plot data obatined from cn_construct_distribution_data
%
% See also cn_construct_distribution_data

% Read input
[ax,args,nargs] = axescheck(varargin{:});
original_args=args;

tint=args{1};
args=args(2:end);

emin=[];

n=1:nargs-1;
for k=n
    specrec{k}=args{k}; 
    % put energy in log eV
    rlog{k} = log10(double(specrec{k}.en))';
    % Pitch angles, turn so that pitch angle 0 is on top
    theta{k} = double(specrec{k}.f)+90;        
    theta_save{k} = double(specrec{k}.f)+90;
    % take out time interval
    if length(tint)==2;
        [~,ind{k}]=irf_tlim(specrec{k}.t,tint);
    elseif length(tint)==1;
        [tt,ind{k}]=irf_tlim(specrec{k}.t,[tint-1 tint+1]);
        ind_min=find(abs((tt-tint))==min(abs((tt-tint))));
        ind={ind{1}(ind_min)};
    end
end

% If no axes is given, initialize figure.
if isempty(ax) 
    ax=irf_plot(1);
end

switch specrec{1}.type
    case 'polar'
        if nargs-1==1; rlog{2}=rlog{1}; end
        % Take out r0
        if isempty(emin) || log10(emin)>min([min(rlog{1}) min(rlog{2})])
            r0 = min([min(rlog{1}) min(rlog{2})]);
        else
            r0 = log10(emin);
        end  

        % Create surf grids    
        for k=n
            r_man{k}=rlog{k}-r0;
            X{k} = r_man{k}*cosd(theta_save{k});
            Y{k} = r_man{k}*sind(theta_save{k});
            C{k} = squeeze(nanmean(log10(specrec{k}.data(ind{k},:,:)),1))';
        end
        if n==1
            r_man{2}=r_man{1};
        end
        if n==1 % Duplicate matrices
            X{2}=X{1};Y{2}=Y{1};C{2}=C{1};
        end            

        if 1 % Plot data
            surf(ax,X{1},Y{1},X{1}*0,C{1}); hold(ax,'on');
            surf(ax,-flipdim(X{2},2),Y{2},X{2}*0,C{2});         
            view(ax,2); axis(ax,'equal','tight'); shading(ax,'flat'); grid(ax,'off');
            cb=colorbar('peer',ax); 
        end
        if 1 % Title and labels
            if length(tint)==2
                t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
                t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
                %prod_str{1}=[t1str,'-',t2str,'UT'];
                %prod_str{2}=product;
                prodstr=[t1str,'-',t2str,' UT    ',specrec{k}.product];
            elseif length(tint)==1
                t1str=datestr(epoch2date(specrec{k}.t(ind{:})),'dd-mmm-yyyy  HH:MM:SS.FFF');
                prodstr=[t1str,' UT    ',specrec{k}.product];
                t2str='';
            end
            %t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            %t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');     
            if n==1 
                str_product=specrec{1}.product;
                str_product(strfind(str_product,'_'))=' ';
                title_str={[t1str,'-',t2str,' UT'],str_product};
                ylabel(cb,specrec{1}.p_label)
            else
                str_product1=specrec{1}.product;
                str_product1(strfind(str_product1,'_'))=' ';
                str_product2=specrec{2}.product;
                str_product2(strfind(str_product2,'_'))=' ';
                title_str={[t1str,'-',t2str,' UT'],...
                    ['Left: ',str_product1],['Right: ', str_product2]};
                ylabel(cb,[char(specrec{1}.p_label), '  /  ', char(specrec{2}.p_label)])
            end
            title(ax,title_str)      
        end
        if 1 % Ticks, something not quite right with these
            xticks=log10([1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]/1000)-r0;
            %x_ind=find(xticks>0);
            xticks=[xticks(find(xticks>0))];
            xticks=[xticks(find(xticks<max([max(r_man{1}) max(r_man{2})])))];
            xticklabels=cell(size(xticks));
            for k=1:length(xticklabels)
                xticklabels{k}=num2str(10.^(xticks(k)+r0)/1000);
            end         
            xticks=[-flipdim(xticks,2) xticks];
            xticklabels=[flipdim(xticklabels,2) xticklabels];
            yticks=xticks;
            yticklabels=xticklabels; 
            set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
                'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)  
            xlabel(ax,'Energy  [keV]'); ylabel(ax,'Energy  [keV]')
        end
        if 1 % Pitch angle labels
            rmax=max([max(r_man{1}) max(r_man{2})]);
            text(0-0.2,rmax-0.5,0,'0^o')
            text(0-0.2,-rmax+0.5,0,'180^o')
            text(-0.2-rmax+0.5,0,0,'90^o')
            text(-0.2+rmax-0.5,0,0,'90^o')
        end            
    case 'cross-section'
        if 1 % Select data to plot
            PA0=squeeze(nanmean(specrec{k}.data(ind{k},1,:),1));
            if mod(length(specrec{k}.f),2)
                PA90=squeeze(nanmean(specrec{k}.data(ind{k},find(length(specrec{k}.f)),:),1));
            else
                PA90=squeeze(nanmean(nanmean(specrec{k}.data(ind{k},[length(specrec{k}.f)/2 length(specrec{k}.f)/2+1],:),1),2));
            end
            PA180=squeeze(nanmean(specrec{k}.data(ind{k},end,:),1));
            PAbg=squeeze(nanmean(nanmean(specrec{k}.bg(ind{k},:,:),1),2)); % One common zero-count level for all levels
        end
        [PA0 PA90 PA180 PAbg]
        if 1 % Plotting data
            loglog(ax,specrec{k}.en',PA0',specrec{k}.en',PA90',specrec{k}.en',PA180',specrec{k}.en',PAbg','--');
            %set(ax,'xlim',[specrec.en(1) specrec.en(end)])
            if length(tint)==2
                t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
                t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
                %prod_str{1}=[t1str,'-',t2str,'UT'];
                %prod_str{2}=product;
                prodstr=[t1str,'-',t2str,' UT    ',specrec{k}.product];
            elseif length(tint)==1
                t1str=datestr(epoch2date(specrec{k}.t(ind{:})),'dd-mmm-yyyy  HH:MM:SS.FFF');
                prodstr=[t1str,' UT    ',specrec{k}.product];
            end
            prodstr(strfind(prodstr,'_'))=' ';
            irf_legend(ax,{'0^o','90^o','180^o','-- Zero count'},[0.94 0.94]) 
            ylabel(ax,specrec{k}.p_label)
            xlabel(ax,'Energy  [eV]')
            grid(ax,'off');
            title(ax,prodstr)
         end
        
end