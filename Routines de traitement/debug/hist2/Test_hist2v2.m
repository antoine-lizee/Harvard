% Test the different results for a range of "softening" values
%% Setting the initial values
clf
f1=2;
f2=2;
sigma=1;
n=4000;
%P=[randn(2000,sigma); simuprob(1000,f1); simuprob(1000,f2)];
P=[randn(2000,sigma); simuprob(2000,f1)];
hold all; 
Ns=cell(1,1);
Xs=cell(1,1);
I=[0,5,150];
LW=[1 2 2];

%% Preparation of the pdf
x=-6:0.1:6;
c1=-cos(x*f1)/2+1/2;
c1=c1/(sum(c1)*12/numel(c1));
c2=-cos(x*f2)/2+1/2;
c2=c2/(sum(c2)*12/numel(c2));
g=exp(-x.^2/(2*sigma^2))*sqrt(1/(2*pi*sigma^2));
% y=c1+c2+g;
% y=y/3;
y=c1+g;
y=y/2;

%% Drawing of the different plots
for i=1:length(I),
    %subplot(1,3,i);
    %figure(i);
    hold all;
    %plot(x,y*n*(range(x)/100),':','DisplayName','Theoric pdf');
    [Ns{i} Xs{i}]=hist2(P,100,I(i));
    plot(Xs{i}, Ns{i}, '-','LineWidth',LW(i),'DisplayName',['softened histogram of level ' num2str(I(i))] );
end

plot(x,y*n*(range(x)/100),':','DisplayName','Theoric pdf');

%% Plot editing
legend(gca, 'show');
