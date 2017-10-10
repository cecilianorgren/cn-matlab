function [dt,dtmax,dtmin]=my_interf(p1,p2,p3,p4)
%
%
% compute time delay between p3 and p4. 
% E3 = p3-(p1+p2)/2
% E4 =-p4+(p1+p2)/2 
E1 = irf_add( 1,p1,-1,irf_add(0.5,p3,0.5,p4));
E2 = irf_add(-1,p2, 1,irf_add(0.5,p3,0.5,p4));
E3 = irf_add( 1,p3,-1,irf_add(0.5,p1,0.5,p2));
E4 = irf_add(-1,p4, 1,irf_add(0.5,p1,0.5,p2));

% This above actually gives electric field looking stuff, although with
% amplitudes of 5000 something.
p3 = irf_tappl(E3,'*0.00212/0.044');
p4 = irf_tappl(E4,'*0.00212/0.044');

p3 = irf_tappl(E3,'/0.044');
p4 = irf_tappl(E4,'/0.044');

MAX1 = 50; %mV/m % CN: all the max should be above this value, but i 
                     % think its a little bit to low so i set it to 60 mV/m
                     % instead, based on:
                     % irf_plot(irf_tappl(E3,'*0.00212/0.044'))
%MAX1 = 60;
MIN1 = 50; % CN: sets MIN1 to something bigger due to overall larger amplitudes
%MIN1 = 40;
DTEFW=0; % p3 and p4 are sampled simultaneously, as there are 2 ADCs.
         % there is a delay between p1/p2 and p3/p4 of 1/36000 sec.
PROBE_DIST = .044; % km

tepoch = p3(:,1);
t = tepoch-tepoch(1); % CN: get time to start at 0.
p3 = p3(:,2); % CN: take away time series.
p4 = p4(:,2);

% find maxima and minima
imap3 = find_max(p3,MAX1); % men vad har 0.25 mV/m med p3 att g?ra, det ?r ju helt olika enheter
imip3 = find_max(-p3,-MIN1);
imap4 = find_max(p4,MAX1);
imip4 = find_max(-p4,-MIN1);
figure;
h(1) = irf_subplot(2,1,-1);
plot(t,p3,'g', t(imap3),p3(imap3),'og', t(imip3),p3(imip3),'dg',...
	t+DTEFW,p4,'b', t(imap4)+DTEFW,p4(imap4),'ob', t(imip4)+DTEFW,p4(imip4),'db')
hold on; grid

h(2) = irf_subplot(2,1,-2);

% find dt
dt = zeros(length(imap3),1);
dtmin = dt; dtmax = dt;
%whos im*
%return;
for i = 1:length(imap3);
    pause
	% p3 is reference, adjust p4 to it
	iminp3 = imip3(imip3>imap3(i)); iminp3 = iminp3(1);
	dAp3 = p3(imap3(i)) - p3(iminp3);
	mp3 = (p3(imap3(i)) + p3(iminp3))/2;
	
	iminp4 = imip4(imip4>imap4(i)); iminp4 = iminp4(1);
	
	% Zero crossings
    % CN: dont know it his works properly, the find_max function is not
    % doing a good job of finding zero-crossnings...
	dAp4 = p4(imap4(i)) - p4(iminp4);
	mp4 = (p4(imap4(i)) + p4(iminp4))/2;
	
	off = fix((iminp4 - imap4(i))/3);
	off = 2;
	datap4 = p4(imap4(i)+off:iminp4-off);
	
	% adjust mean for
    
    p4
	datap4 = datap4 - mp4;
	datap3 = p3(imap3(i)+off:iminp3-off) -mp3;

	% rescale p4
	datap4 = datap4*dAp3/dAp4;
	plot(t(imap3(i)+off:iminp3-off),datap3,'g',...
		t(imap4(i)+off:iminp4-off) +DTEFW, datap4,'b')
	hold on
	pp3 = polyfit(t(imap3(i)+off:iminp3-off), datap3,1);
	pp4 = polyfit(t(imap4(i)+off:iminp4-off), datap4,1);
	plot(t(imap4(i)+off:iminp4-off)+DTEFW,pp4(1)*t(imap4(i)+off:iminp4-off)+pp4(2),'b--',...
		t(imap3(i)+off:iminp3-off), pp3(1)*t(imap3(i)+off:iminp3-off)+pp3(2),'g--')
	tp40 = -pp4(2)/pp4(1); tp30 = -pp3(2)/pp3(1);
	dt(i) = tp40 + DTEFW - tp30 ; % p4 is sampled DTEFW sec later then p3
	plot(tp30,0,'g+',tp40+DTEFW,0,'b+')
	axes(h(1))
	text(tp30,0,sprintf('%.2fe3',PROBE_DIST/1000/dt(i)))
	fprintf('%s    dt=%f ms v=%.2fe3 km/s \n',epoch2iso(tepoch(1)+tp30,1),dt(i)*1000,PROBE_DIST/1000/dt(i))
	axes(h(2))

	% Maxima
	off = 2;
	
	plot(t(imap3(i)-off:imap3(i)+off), p3(imap3(i)-off:imap3(i)+off),'g',...
		t(imap4(i)-off:imap4(i)+off)+DTEFW, p4(imap4(i)-off:imap4(i)+off),'b')
	
	pp3 = polyfit(t(imap3(i)-off:imap3(i)+off), p3(imap3(i)-off:imap3(i)+off),2);
	pp4 = polyfit(t(imap4(i)-off:imap4(i)+off), p4(imap4(i)-off:imap4(i)+off),2);
	
	plot(t(imap4(i)-off:imap4(i)+off)+DTEFW,...
		pp4(1)*t(imap4(i)-off:imap4(i)+off).^2+pp4(2)*t(imap4(i)-off:imap4(i)+off)+pp4(3),'b--',...
		t(imap3(i)-off:imap3(i)+off),...
		pp3(1)*t(imap3(i)-off:imap3(i)+off).^2+pp3(2)*t(imap3(i)-off:imap3(i)+off)+pp3(3),'g--')
	
	tp40 = -pp4(2)/pp4(1)/2; tp30 = -pp3(2)/pp3(1)/2;
	
	dtmax(i) = tp40 + DTEFW - tp30 ; % p4 is sampled 1/36000 sec later then p3
	plot(tp30,p3(imap3(i)),'g+',tp40+DTEFW,p4(imap4(i)),'b+')
	axes(h(1))
	text(tp30,p3(imap3(i)),sprintf('%.2fe3',PROBE_DIST/1000.0/dtmax(i)))
	disp(sprintf('%s dtmax=%f ms v=%.2fe3 km/s',epoch2iso(tepoch(1)+tp30,1),dtmax(i)*1000.0,PROBE_DIST/1000.0/dtmax(i)))
	axes(h(2))
	
	% Minima
	plot(t(iminp3-off:iminp3+off), p3(iminp3-off:iminp3+off),'g',...
		t(iminp4-off:iminp4+off)+DTEFW, p4(iminp4-off:iminp4+off),'b')
	
	pp3 = polyfit(t(iminp3-off:iminp3+off), p3(iminp3-off:iminp3+off),2);
	pp4 = polyfit(t(iminp4-off:iminp4+off), p4(iminp4-off:iminp4+off),2);
	
	plot(t(iminp4-off:iminp4+off)+DTEFW,...
		pp4(1)*t(iminp4-off:iminp4+off).^2+pp4(2)*t(iminp4-off:iminp4+off)+pp4(3),'b--',...
		t(iminp3-off:iminp3+off),...
		pp3(1)*t(iminp3-off:iminp3+off).^2+pp3(2)*t(iminp3-off:iminp3+off)+pp3(3),'g--')
	
	tp40 = -pp4(2)/pp4(1)/2; tp30 = -pp3(2)/pp3(1)/2;
	
	dtmin(i) = tp40 + DTEFW - tp30 ; % p4 is sampled 1/36000 sec later then p3
	plot(tp30,p3(iminp3),'g+',tp40+DTEFW,p4(iminp4),'b+')
	axes(h(1))
	text(tp30,p3(iminp3),sprintf('%.2fe3',PROBE_DIST/1000/dtmin(i)))
	disp(sprintf('%s dtmin=%f ms v=%.2fe3 km/s',epoch2iso(tepoch(1)+tp30,1),dtmin(i)*1000,PROBE_DIST/1000/dtmin(i)))
	axes(h(2))
end	
hold off; grid

function imaxima = find_max(sig,maxv)
% sig is like -8000 while maxv is 0.025, so condition on third line always
% holds and all is empty, so it doesnt work.
imaxima = [];
ndata = length(sig);
sig(sig<maxv) = NaN;
i1 = 1;
iter=0;
while 1
	i1_tmp = find(~isnan(sig(i1:end))); % find all indices that is above maxv
	if isempty(i1_tmp), break, end
	i1 = i1 + i1_tmp(1);
	i2 = find(isnan(sig(i1+1:end)))-1;
	if isempty(i2), i2 = ndata;
	else i2 = i2(1) + i1;
	end
	im = find(sig(i1:i2) == max(sig(i1:i2)));
    if numel(im)>1 % takes first index if there are two with exact values
        im=im(1);
    end
	imaxima = [imaxima i1+im-1];
	i1 = i2 + 1;
	if i1>=ndata, break, end
    %iter=iter+1;
    %disp(num2str(iter))
    %if 0 ; %iter==228
    %    disp('Here comes problem. Debug.')
    %end
end