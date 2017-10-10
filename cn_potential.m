function phi= cn_potential(E,ecoord,v,vcoord)
ecols=size(E,2);
erows=size(E,1);
vrows=size(v,1);
vcols=size(v,2);

if vrows == 1
    switch vcols
        case 3
            v=c_coord_trans(vcoord,ecoord,[E(1,1) v],'cl_id',3);
            vvec=ones(erows,4);
            vvec(:,2)=v(2);
            vvec(:,3)=v(3);
            vvec(:,4)=v(4);
        case 4        
            v=c_coord_trans(vcoord,ecoord,v,'cl_id',3);
            vvec(:,2)=v(2);
            vvec(:,3)=v(3);
            vvec(:,4)=v(4);
    end
else
            vvec=c_coord_trans(vcoord,ecoord,v,'cl_id',3);
            vvec(:,1)=1;
end

phit=irf_integrate(E,E(end,1));      
phi=phit.*vvec(:,1:ecols);
phi=[phi(:,1) -sum(phi(:,2:end),2)];