model = 1;
figure(12);
set(gcf,'position',[206 1 565 923]);
set(gcf,'position',[206 1 923 565]);
for model=1:2;         
    switch model
        case 1 % 2007-08-31
            lr = 15;
            lz = 10;
            phi0 = 200;
            B0 = 25;
            Tpar = 1600;
            Tper = 1600;
            r0 = 15;
            a1 = 45;
            a2 = 0;        
        case 2 % Tao2011
            lr = 35;
            lz = 35;
            phi0 = 3000;
            B0 = 50;
            Tpar = 8000;
            Tper = 5000;
            r0 = lr;
            a1 = 45;
            a2 = 0;                       
    end

    subplot(2,1,model)
    [h,x,y,z,vx,vy,vz,ax,ay,az]=art2.ExB_plot(lr,lz,phi0,B0,Tpar,Tper,r0,a1,a2);
    axel(model) = h;
end