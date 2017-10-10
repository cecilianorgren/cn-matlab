function [side, area, mperpixel] = uranus_pic_size(instrument, distance)
Ur=25559; % km
mperpixel=0;
switch instrument
    case 'CAM'       
        res=300; %m/pixel
        pixels=1024;
        deg=0.29;
        side=distance*atand(deg);
        mperpixel=side/pixels;
        area=side.*side;
    case 'VIR'
        res=0;
        mrad=64e-3;
        side=distance*atan(mrad);
        area=side.*side;
    case 'TIR-far'
        mrad=4e-3;
        side=distance*atan(mrad);
        area=2*pi*side;
    case 'TIR-mid'
        mrad=0.3e-3;
        side=distance*atan(mrad);
        area=2*pi*side;
end
