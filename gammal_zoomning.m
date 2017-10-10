    switch zoom % zoom 20070831 
        case 1
           irf_zoom(h,'x',tint(9,:)); 
        case 2
           irf_zoom(h,'x',[toepoch([2007 09 02 14 30 30]) toepoch([2007 09 02 14 36 00])]);
        case 3
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 00]) toepoch([2007 09 02 14 33 00])]);
        case 4
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 20]) toepoch([2007 09 02 14 32 45])]);
        case 5
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 26]) toepoch([2007 09 02 14 32 32])]);
    end
    switch zoom % zoom 20070902 
        case 1
           irf_zoom(h,'x',tint(10,:)); 
        case 2
           irf_zoom(h,'x',[toepoch([2007 09 02 14 30 30]) toepoch([2007 09 02 14 36 00])]);
        case 3
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 00]) toepoch([2007 09 02 14 33 00])]);
        case 4
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 20]) toepoch([2007 09 02 14 32 45])]);
        case 5
           irf_zoom(h,'x',[toepoch([2007 09 02 14 32 26]) toepoch([2007 09 02 14 32 32])]);        
    end   
    switch zoom % zoom 20070926 