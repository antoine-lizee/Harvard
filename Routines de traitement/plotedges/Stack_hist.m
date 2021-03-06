% Plot all the microsorum and the Zinnia in the same Histograms :

figure(20); clf; h1=newplot;
hold all;
figure(21); clf; h2=newplot;
hold on;
linspec={':' ':' ':' ':' ':' '-'};
LW=[1 1 1 1 1 2];


for i=1:6    
    Ntot=length(p(i).t);
    hist_graph=p(i).t(:)*180/pi;
    hist_graph(hist_graph==0)=[];
    
    [Na(i,:) xouta]=hist2(hist_graph,100,2);
    Na(i,:)=Na(i,:)*100/Ntot; % Convert in "%" for each edge (bottom and top)
    plot(h1,xouta, Na(i,:),'-','LineWidth',LW(i));
    title(h1,{'Histogram of the angles of the up and down edges';'Comparison between microsorum and zinnia'});
    
    [Ns xouts]=hist2(p(i).r,20);
    Ns=Ns*100/Ntot; % Convert in "%" for each edge (bottom, top and right)
    plot(h2,xouts, Ns,linspec{i},'MarkerSize',4);
%     legend('bottom edge', 'right edge', 'top edge');
    title(h2,{'Histogram of the edges size';'Comparison between microsorum and zinnia'});

end