n=1000;

figure(204); clf; 
h4=newplot;
hold all;

sigma=1;
P=randn(n,sigma);

N=cell(5,1);
X=cell(5,1);

[N{6} X{6}]=histv2(P,10);
plot(X{6},N{6},'DisplayName','normal 10-bin histogram','Linewidth',1);
%bar(X{6},N{6},'DisplayName','normal 10-bin histogram','Linewidth',1);

%[N{1} X{1}]=histc(P,e);
%plot(X{1},N{1}*2);

% [N{2} X{2}]=histv3(P,10,50);
% plot(X{2},N{2});
% 
% [N{3} X{3}]=histv4(P,10,50);
% plot(X{3},N{3},'Linewidth',0.3);

[N{4} X{4}]=hist(P,50);
plot(X{4},N{4}*5,'DisplayName','normal 50-bin histogram');

[N{5} X{5}]=histv2(P,10,50);
h=plot(X{5},N{5},'DisplayName',sprintf('window histogram with \n 10th-wide bins but 50 points'),'Linewidth',2);
uistack(h,'bottom');

x=-4:0.1:4;
y=exp(-x.^2/(2*sigma^2))*sqrt(1/(2*pi*sigma^2));
plot(x,y*n*(range(P)/10),':','DisplayName','gaussian pdf');

title({['Comparison between classic histograms and window histogram : Gaussian exemple (' num2str(n) 'points)']});

legend1 = legend(gca,'show');

