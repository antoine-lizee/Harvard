clf
    hold on; 
    [Ns2 xouts2]=histv4(p.r,35,70);
    plot(xouts2, Ns2, '-');
    
    [Ns xouts]=hist2(p.r,60,4);
    plot(xouts, Ns, '-');
    
    [Ns xouts]=hist2(p.r,60,5);
    plot(xouts, Ns, '-');
    
     [Ns xouts]=hist2(p.r,60,6);
    plot(xouts, Ns, '-');
    
    [Ns xouts]=hist2(p.r,60,3);
    plot(xouts, Ns, '-');