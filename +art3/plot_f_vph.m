% plots the electron distributions and mark out the phase velocities

e = 1.6022e-19;
kB = 1.38e-23;
me = 9.1094e-31;
m = 9.1094e-31;
mp = 1.6726e-27;

% f_e
f_int = @(v,vt,n,vd) n./(sqrt(pi)*vt).*exp(-(v-vd).^2/vt/vt);
f_3D = @(v,vt,n,vd) n./(sqrt(pi)*vt).^3.*exp(-(v-vd).^2/vt/vt);

f = f_int;
n = 1*1e6;
Tebg = 300;
Tebeam = 12;
vtbg = sqrt(2*qe*Tebg/me);
vtbeam = sqrt(2*qe*Tebeam/me);
v = linspace(-2*vtbg,2*vtbg,1000);
R=s{4}.R;
nR = numel(R);
nS = numel(S);
for iR=11:nR;
    for iS=1:nS;
        subplot(2,1,1)
        f=f_int;
        semilogy(v/vtbg,f(v,vtbg,n*(1-R(iR)),0)+f(v,vtbeam,n*R(iR),S(iS)*vtbg));
        vTvtbeam=(S(iS)*vtbg-s{4}.vph(iR,iS))/vtbeam;
        title(['R=' num2str(R(iR)) ', S=' num2str(S(iS)) ', v_T/v_{t,beam} = ' num2str(vTvtbeam,'%.1f')])
        %hold on; plot(vtbg*[1 1],get(gca,'ylim'),'g--'); hold off
        hold on; plot(S(iS)*[1 1]-vtbeam/vtbg,get(gca,'ylim'),'g--'); hold off;
        try
            hold on;
            plot(s{4}.vph(iR,iS)/vtbg*[1 1],get(gca,'ylim'),'r--')
            hold off
        end
        hold on; hl(1) = plot(3*vtbg,mean(get(gca,'ylim')),'r--'); hl(2) = plot(3*vtbg,mean(get(gca,'ylim')),'g--'); hold off
        legend([hl(1:2)],'v_{ph}','v_{beam}-v_{te,beam}','location','eastoutside')
        set(gca,'xlim',2*[-1 1])
        xlabel('v/v_{t,bg}'); ylabel('\int f(v_x,v_y,v_z)dydz')
        
        subplot(2,1,2)
        f=f_3D;
        semilogy(v/vtbg,f(v,vtbg,n*(1-R(iR)),0)+f(v,vtbeam,n*R(iR),S(iS)*vtbg));
        title(['R=' num2str(R(iR)) ', S=' num2str(S(iS)) ', v_T/v_{t,beam} = ' num2str(vTvtbeam,'%.1f')])
        %hold on; plot(vtbg*[1 1],get(gca,'ylim'),'g--'); hold off
        hold on; plot(S(iS)*[1 1]-vtbeam/vtbg,get(gca,'ylim'),'g--'); hold off;
        try
            hold on;
            plot(s{4}.vph(iR,iS)/vtbg*[1 1],get(gca,'ylim'),'r--')
            hold off
        end
        hold on; hl(1) = plot(3*vtbg,mean(get(gca,'ylim')),'r--'); hl(2) = plot(3*vtbg,mean(get(gca,'ylim')),'g--'); hold off
        legend([hl(1:2)],'v_{ph}','v_{beam}-v_{te,beam}','location','eastoutside')        
        set(gca,'xlim',2*[-1 1])
        xlabel('v/v_{t,bg}'); ylabel('f(v_x,0,0)')
        
        cn.print(['f_vph_iR' num2str(iR) 'iS' num2str(iS)])
        close
    end
end
    