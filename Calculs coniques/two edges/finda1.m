findx = @(x,y)...
cos(x)-cos(pi/2-y-x)+(sin(x)*(1/tan(y)+1/(sqrt(2)*sin(y)))-sin(pi/2-x-y)*(1/tan(y)+sqrt(2)/sin(y)));

figure(200);
clf;
hold all;
for alpha0=0:5:180
    alpha=pi/180*alpha0+2*eps;
    a1=0:0.01:(90-alpha0)*pi/180;
    f=findx(a1,alpha);
    plot(a1,f,'DisplayName',['alpha0=' num2str(alpha0)]);
end
Xtick=(0:30:91)*pi/180;
set(gca,'XTick',Xtick)
set(gca,'XTickLabel',{'0','30','60','90'})
legend1=legend('show');
set(legend1,'Location','NorthEastOutside');
plot([0,pi/2],[0, 0]);

alpha0=0:0.01:180;
a1=zeros(1,length(alpha0));
alpha=pi/180*alpha0+2*eps;
for i=1:length(alpha0);
    options=optimset('Display', 'notify');
    [a1(i) fval exitflag]=fzero(@(x) findx(x,alpha(i)), 0.5, options);
end
a1=a1*180/pi;
figure(201);
plot(alpha0,a1);
results=[alpha0;a1];
disp('pour un angle de 65° :');
disp(results(2,results(1,:)==65));
results_geogebra=results(:,1:10:end)';
