txtfile  = '/Users/cecilia/Data/PIC/data/fields-01200.dat';

[xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
    jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz, ...
    dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a, ...
    wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez] ...
    = func_read_file_fields_oxygen(txtfile);