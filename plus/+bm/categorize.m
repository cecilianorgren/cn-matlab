% Determine duration and location of burst modes
% Run bm.list before to obtain BM1s and BM2s

% BM1
if 0
    tintISO = regexp(BM1s,'(\d*)-(\d*)-(\d*)T(\d*):(\d*):(\d*)Z','match');
    nBM1 = numel(BM1s);
    for k = 1:nBM1
        tint(k,:) = irf_time([tintISO{k}{1}; tintISO{k}{2}],'iso2epoch')';
    end    
else
    bm1TT = irf.TimeTable('/Users/Cecilia/Data/BM1.txt');
    tint = bm1TT.TimeInterval;    
end
T = tint(:,2) - tint(:,1);

if 0 % plot time duration of BM
    subplot(2,1,1)
    hist(T/60/60,20)
    xlabel('T_{BM} [h]'); ylabel('Occurence')
    subplot(2,1,2)
    plot(tint(:,1),T/60/60,'o')
    xlabel('Time [epoch]'); ylabel('T_{BM} [h]')
end

% Select the events that are in the magnetotail box
c_eval('tailC?TT=irf.tt(''read_IRF'',''C?_in_tailbox'');',1:4);
% Events that are in BM and magnetotail box, one for each SC
c_eval('tailBMC? = intersect(tailC?TT,bm1TT);',1:4);
% Any SC being in tailbox
% write new class function for this...

%% List dates when given satellite was in tailbox
sc = 4;
c_eval('nEvent = numel(tailBMC?);',sc);
c_eval('irf_time(tailBMC?.TimeInterval,''yyyymmdd'')',sc);

