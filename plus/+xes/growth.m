function growthrate = growth(t1,E1,t2,E2)

growthrate = log(E1./E2)./(t1-t2);