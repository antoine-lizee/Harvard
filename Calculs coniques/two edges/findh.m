findh1 = @(H1,H2,alpha)...
(H2./sin(alpha)-H1./tan(alpha));

figure(201);
clf;
hold all;
H2=2;
for H1=1:0.2:3
    alpha=0:0.01:90*pi/180;
    h1=findh1(H1,H2,alpha);
    plot(alpha,h1,'DisplayName',['H1=' num2str(H1)]);
end
Xtick=(0:15:91)*pi/180;
set(gca,'XTick',Xtick)
set(gca,'XTickLabel',{'0','15','30','45','60','75','90'})
legend1=legend('show');
set(legend1,'Location','NorthEastOutside');
axis([-Inf Inf -3 4])