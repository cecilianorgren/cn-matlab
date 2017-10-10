function phipar = fitparpotential(ne,Te,PSD0,PSD180,energy)

theta = [7.5 172.5];
Evec = energy;

me = 9.1e-31;
qe = 1.6e-19;

SHne = ne;
SHTe = Te;

phipar = [0:1:200];
phiparfit = zeros(length(phipar),1);

dist1fit = SHne*(me/(2*pi*qe*SHTe))^(3/2)*exp(-Evec/SHTe);

trapdist = zeros(length(phipar),length(Evec),length(theta));

for ii=1:length(phipar)
    for jj=1:2
        for kk=1:length(Evec)
            Epar = Evec(kk)*cosd(theta(jj))^2;
            if(Epar >= phipar(ii)) 
                trapdist(ii,kk,jj) = SHne*(me/(2*pi*qe*SHTe))^(3/2)*exp(-(Evec(kk) - phipar(ii))/SHTe);
            else 
                trapdist(ii,kk,jj) = SHne*(me/(2*pi*qe*SHTe))^(3/2)*exp(-(Evec(kk)*sind(theta(jj))^2)/SHTe);
            end
        end
    end
end

for ii=1:length(phipar)
    phiparfit(ii) = sum(abs(PSD0-log10(trapdist(ii,:,1))).^2)+sum(abs(PSD180-log10(trapdist(ii,:,2))).^2);
end


[~,minpos] = min(phiparfit);
phipar = phipar(minpos);

%loglog(energy,PSD0,energy,PSD180)
end