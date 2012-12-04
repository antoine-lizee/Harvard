%If launched as opened, this script will compare different results for
%increasing steps of softening and decreasing values of initial bin width.

%linespec={'-r','-g', '-b','-y', '-k', '-w'};
figure(201); clf; 
h1=newplot;
hold all;
figure(202); clf; 
h2=newplot;
hold all;
figure(203); clf; 
h3=newplot;
hold all;

Nst=cell(2,2);
Xout=cell(2,2);

%P=p(6).r;
P=randn(1,1000);
i1=[1,4,10,20,50]
i2=1:2;
for i=1:length(i1)
    [Nst{1,i} Xout{1,i}]=hist2(P,20+40*(i1(i)-1),i1(i)-1);
    Nst{1,i}=Nst{1,i}*(1+2*(i1(i)-1));
    plot(h3,Xout{1,i},Nst{1,i});
end
legend(h3,num2str(i1'))

test=[0.501,0.01,0.99];
c={'r','g','b','y','c','m','k'};
for i=1:length(i1)
    [Nst{3,i} Xout{3,i}]=hist2(test,20+40*(i1(i)-1),i1(i)-1);
    Nst{3,i}=Nst{3,i}*(1+2*(i1(i)-1));
    %bar(h1,Xout{3,i},Nst{3,i},'LineStyle', 'none','FaceColor',c{i},'BarWidth',0.5);
    plot(h1,Xout{3,i},Nst{3,i});
end
legend(h1,num2str(i1'))

for i=1:length(i2)
    [Nst{2,i} Xout{2,i}]=hist(P,20+40*(i2(i)-1));
    Nst{2,i}=Nst{2,i}*(1+2*(i2(i)-1));
    plot(h2,Xout{2,i},Nst{2,i});
end
legend(h2,num2str(i2'))