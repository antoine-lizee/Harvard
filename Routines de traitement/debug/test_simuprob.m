N=10000;
f=4;
[ P ] = simuprob( N, f);
hold all
hist(P,100);

x=-6:0.1:6;
c1=-cos(x*f1)/2+1/2;
plot(x,c1*N*(12/100)/(sum(c1)*12/numel(c1)),'Linewidth',3);